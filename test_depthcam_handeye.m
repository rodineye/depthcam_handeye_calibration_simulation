% the true value of hand-eye calibration
% [Rc Tc; 0 1]
Tc = [0.0325;0.27;0.0438];
Rc = [1 0 0; 0 -1 0; 0 0 -1];
%Rc = [0.9942    0.0537   -0.0931;
%      0.0398   -0.9886   -0.1450;
%     -0.0998    0.1404   -0.9850];
invRc = inv(Rc); 

% the true value of TCP calibration
TCP = [0;0;0.345];

%%%%%
data_n = 30;
count = 50;
sum_r = [0 0 0];
sum_tc = [0; 0; 0];

for ni = 1:count
    [Rc_1 Tc_1] = depthcam_handeye(data_n, TCP, Rc, Tc);
    [r1 r2 r3] = dcm2angle(invRc*Rc_1,'xyz');
    all_r(ni) = norm([r1 r2 r3]);
    all_tc(ni) = norm(Tc_1 - Tc);
    sum_r = sum_r + [r1 r2 r3];
    sum_tc = sum_tc + Tc_1;
end

sum_r = sum_r / count;

sum_tc = sum_tc / count - Tc;
   
figure;
plot((1:count), all_r, 'k-', (1:count), all_tc, 'k-.');
title('50 measurement error distribution');
xlabel('Nth measurement');
ylabel('error');
legend('orientation error','position error');

Rc_1 = Rc*angle2dcm(sum_r(1),sum_r(2),sum_r(3),'xyz')
       
Tc_1 = sum_tc + Tc
       
%%%%%%%%%%%%%
clear all;

% the true value of hand-eye calibration
% [Rc Tc; 0 1]
Tc = [0.0325;0.27;0.0438];
Rc = [1 0 0; 0 -1 0; 0 0 -1];
%Rc = [0.9942    0.0537   -0.0931;
%      0.0398   -0.9886   -0.1450;
%     -0.0998    0.1404   -0.9850];
invRc = inv(Rc); 

% the true value of TCP calibration
TCP = [0;0;0.345];


% Set the number of TCP contact points and hand-eye sampling poses to 10, 15, 20, 25, 
% and 30 respectively, and perform simulation calculations.  
data_n = [10 15 20 25 30];
count = 50;
for nj = 1:5
    sum_r = [0 0 0];
    sum_tc = [0; 0; 0];

    for ni = 1:count
        [Rc_1 Tc_1] = depthcam_handeye(data_n(nj), TCP, Rc, Tc);
        [r1 r2 r3] = dcm2angle(invRc*Rc_1,'xyz');
        sum_r = sum_r + [r1 r2 r3];
        sum_tc = sum_tc + Tc_1;
    end

    sum_r = sum_r / count;

    sum_tc = sum_tc / count - Tc;

    All_ori1(nj) = norm(sum_r);

    All_dis1(nj) = norm(sum_tc);
end

figure;
plot(data_n, All_ori1, 'k-', data_n, All_dis1, 'k-.');
title('The influence of the number of samples on the calibration accuracy');
xlabel('the number of samples');
ylabel('error');
legend('orientation error','position error');


