

% "High precise and zero-cost solution for fully automatic industrial robot TCP calibration" 

function [T] = tcp_calib(data,n)

    %plane norm vector
    normal_vector = [n(1) n(2) -1];
    %normalize plane norm vector
    nor = norm(normal_vector);
    normal_vector = [n(1)/nor n(2)/nor -1/nor];
   
    %paper algorithm 2
    len = length(data)/3; 
    M = [];
    D = [];
    for i = 1:len
        M(i,1:3) = normal_vector*data(3*i-2:3*i,1:3);
        D(i,1) = -normal_vector*data(3*i-2:3*i,4);
    end
    k = 1;
    A = [];
    B = [];
    for i = 1:(len-1)
        for j = i+1:len
            A(k,1:3) = M(i,1:3)-M(j,1:3);
            B(k,1) = D(i,1)-D(j,1);
            k = k+1;
        end
    end

    [Q,R] = qr(A);
    T = R\(Q'*B);
   
    %T = A \ B
   
end