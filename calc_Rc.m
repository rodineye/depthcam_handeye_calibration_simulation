
% calc the reprojection error to obtain Rc

function Rc = calc_Rc(A,B)

    min_sum_dis = 1000000;
    
    for ni = 1:length(A)
        for nj = ni+1:length(A)
            for nk = nj+1:length(A)
               
                A1(1:3,1) = A(1:3,ni);
                A1(1:3,2) = A(1:3,nj);
                A1(1:3,3) = A(1:3,nk);
               
                B1(1:3,1) = B(1:3,ni);
                B1(1:3,2) = B(1:3,nj);
                B1(1:3,3) = B(1:3,nk);
               
                % formula 12
                Rc_1 = B1*inv(A1);
                
                % formula 13               
                cur_sum_dis = sum(sum(abs(Rc_1*A-B)));               
                if (cur_sum_dis < min_sum_dis)
                    min_sum_dis = cur_sum_dis;
                    Rc = Rc_1;
                end
               
            end
        end
    end

end