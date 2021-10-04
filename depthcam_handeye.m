function [calib_Rc calib_Tc] = depthcam_handeye(data_n, TCP, Rc, Tc)

% In the simulation, the error of the depth value (z direction) obtained by the depth camera is set to be a random error between -3mm and 3mm, 
% and the xy direction is a random error between -1.5mm and 1.5mm (in accordance with the performance of Intel realsense d435 camera used in section 4.2 
% when the object distance is 25cm~50cm)
x_noise_level = 0.0015;
y_noise_level = 0.0015;
z_noise_level = 0.003;

% import calibration plane points
% the measurement error of position of contact point is 0~0.2mm
%data = load('normalsim0201.txt');
% the measurement error of position of contact point is 0~0.5mm
data = load('normalsim0501.txt');

% data_true is true data
data_true = data(1:data_n,1:3);
% data_noisy is noisy data
data_noisy = data(1:data_n,4:6);

len = data_n;

% fit plane with true data
true_plane_equ = plane_fit(data_true);
true_plane_norm_vec = [true_plane_equ(1),true_plane_equ(2),-1];
true_plane_equ(3) = true_plane_equ(3) / norm(true_plane_norm_vec);
true_plane_norm_vec = true_plane_norm_vec / norm(true_plane_norm_vec);

% fit plane with noisy data
noisy_plane_equ = plane_fit(data_noisy);
noisy_plane_norm_vec = [noisy_plane_equ(1),noisy_plane_equ(2),-1];
noisy_plane_equ(3) = noisy_plane_equ(3) / norm(noisy_plane_norm_vec);
noisy_plane_norm_vec = noisy_plane_norm_vec / norm(noisy_plane_norm_vec);

% import Rb
Homo = load('homo.txt');
Homo = Homo(1:data_n*3,:);

% calc Tb and updata Homo
point = data_true';
for i = 1:len
    Rb = Homo(i*3-2:i*3,1:3);
    Homo(i*3-2:i*3,4) = point(:,i) - Rb * TCP;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TCP calibration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute noisy Tb
point = data_noisy';
for i = 1:len
    Rb = Homo(i*3-2:i*3,1:3);
    Tb = point(:,i) - Rb * TCP;
    noisy_Homo(i*3-2:i*3,1:3) = Rb;
    noisy_Homo(i*3-2:i*3,4) = Tb;
end

% TCP calibration
TCP_noisy = tcp_calib(noisy_Homo, noisy_plane_equ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%hand-eye calibration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obtain camera sim data
for j = 1:len
    Rb = Homo(j*3-2:j*3,1:3);
    Tb = Homo(j*3-2:j*3,4);
   
    % compute To
    % formula 17, AX = B
    % Rb*Rc*To + Rb*Tc + Tb = [x y z]'
    for i = 1:len
        % true value
        To = (Rb*Rc) \ (data_true(i, 1:3)' -Rb*Tc -Tb);
        % add noise
        To_noise(i,1:3) = To' + [-x_noise_level+x_noise_level*2*rand() -y_noise_level+y_noise_level*2*rand() -z_noise_level+z_noise_level*2*rand()];
        All_To_noise(j,i,1:3) = To_noise(i,1:3);
    end
end

% frist phase, compute Rc
for j = 1:len
    Rb = noisy_Homo(j*3-2:j*3,1:3);
   
    for i = 1:len
        To_noise(i,1:3) = All_To_noise(j,i,1:3);
    end
   
    % camera plane a'x + b'y - z + d' = 0
    cam_plane_equ = plane_fit(To_noise);
    cam_plane_norm_vec = [cam_plane_equ(1),cam_plane_equ(2),-1];
    cam_plane_norm_vec = cam_plane_norm_vec / norm(cam_plane_norm_vec);

    % Rb*Rc*cam_plane_norm_vec' = noisy_plane_norm_vec'
    % dcm*cam_plane_norm_vec' = noisy_plane_norm_vec'
    % formula 8,9
    rot_axis = cross(cam_plane_norm_vec,noisy_plane_norm_vec);
    rot_quat = [1 + dot(cam_plane_norm_vec,noisy_plane_norm_vec),rot_axis(1),rot_axis(2),rot_axis(3)];
    rot_quat = quatnormalize(rot_quat);   

    % quat2dcm is different from matlab
    qin = rot_quat;
    dcm = zeros(3,3,size(qin,1));
    dcm(1,1,:) = qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2;
    dcm(1,2,:) = 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4));
    dcm(1,3,:) = 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3));
    dcm(2,1,:) = 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4));
    dcm(2,2,:) = qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2;
    dcm(2,3,:) = 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2));
    dcm(3,1,:) = 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3));
    dcm(3,2,:) = 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2));
    dcm(3,3,:) = qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2;

    % verify inv(Rb)*dcm*cam_plane_norm_vec' = Rc*cam_plane_norm_vec'
    % formula 10,11
    v_left = inv(Rb)*dcm*cam_plane_norm_vec';
    %v_right = Rc*cam_plane_norm_vec';
   
    A1(1:3,j) = cam_plane_norm_vec';
    B1(1:3,j) = v_left(1:3,1);
end

% calc the reprojection error to obtain Rc
calib_Rc = calc_Rc(A1,B1);

% formula 14,15
[U,S,V] = svd(calib_Rc); 
calib_Rc = sign(det(S))*U*V';

% second phase, compute Tc
for j = 1:len
    Rb = noisy_Homo(j*3-2:j*3,1:3);
    Tb = noisy_Homo(j*3-2:j*3,4);
   
    for i = 1:len
        To_noise(i,1:3) = All_To_noise(j,i,1:3);
    end

    % formula 19
    % AX = B       
    for i = 1:len
        A2(i,1:3) = noisy_plane_norm_vec*Rb;
        B2(i,1) = -noisy_plane_equ(3) - noisy_plane_norm_vec*Tb - noisy_plane_norm_vec*Rb*calib_Rc*To_noise(i,1:3)';
    end
   
    total_A2(j*len-len+1:j*len,1:3) = A2;
    total_B2(j*len-len+1:j*len,1) = B2;     
end

% formula 20
calib_Tc = total_A2 \ total_B2;

end
