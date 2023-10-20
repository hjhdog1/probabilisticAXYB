function H12 = buildH12(A,B,alpha)
% A,B : SE(3)
    N = size(A,3);
    H12 = zeros(18,6);
    
    for i = 1:N
        b = B(1:3,4,i);
        b_hat = zeros(9,3);
        b_hat(1:3,1) = b;
        b_hat(4:6,2) = b;
        b_hat(7:9,3) = b;
        
        H12(10:18,1:3) = H12(10:18,1:3) - b_hat*A(1:3,1:3,i);
        H12(10:18,4:6) = H12(10:18,4:6) + b_hat;
    end
    
    H12 = alpha*H12;

end