function f = buildf(A,B,alpha)
% A,B : SE(3)
    N = size(A,3);
    f = zeros(24,1);
    
    for i = 1:N        
        b = B(1:3,4,i);
        b_hat = zeros(9,3);
        b_hat(1:3,1) = b;
        b_hat(4:6,2) = b;
        b_hat(7:9,3) = b;
                
        f(10:18,1) = f(10:18,1) - b_hat*A(1:3,4,i);
        f(19:21,1) = f(19:21,1) + A(1:3,1:3,i)'*A(1:3,4,i);
        f(22:24,1) = f(22:24,1) - A(1:3,4,i);
    end
    f = alpha*f;

end