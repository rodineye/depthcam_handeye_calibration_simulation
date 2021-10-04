
% apply the least square method to calc the plane equation
% formula 2

function norm_vector = plane_fit(data)

    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    n = length(data); 

    x_sum = sum(x);
    y_sum = sum(y);
    z_sum = sum(z);
    xx_sum = sum(x.*x);
    yy_sum = sum(y.*y);
    zz_sum = sum(z.*z);
    xy_sum = sum(x.*y);
    xz_sum = sum(x.*z);
    yz_sum = sum(y.*z);

    % ax + by - z + d = 0  
    % AX = B
    A = [xx_sum xy_sum x_sum;
        xy_sum yy_sum y_sum;
        x_sum y_sum n];
    B = [xz_sum;yz_sum;z_sum];

    coff = A\B ;
    a = coff(1);
    b = coff(2);
    d = coff(3);

    % verify plane equation
    %test = a*x + b*y - z + d;

    %plane norm vector [a b -1];
    norm_vector = [a b d];%平面法向量归一化

end