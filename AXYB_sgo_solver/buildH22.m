function H22 = buildH22(A,alpha)
% A,B : SE(3)
    N = size(A,3);
    H22 = N*eye(6,6);
    
    for i = 1:N
        H22(4:6,1:3) = H22(4:6,1:3) - A(1:3,1:3,i);
    end
    H22(1:3,4:6) = H22(4:6,1:3)';
    H22 = alpha*H22;

end