function r = buildr(A,B,alpha)
% A,B : SE(3)
    N = size(A,3);
    r = 0;

    for i = 1:N
        a = A(1:3,4,i);
        b = B(1:3,4,i);
        r = r + a'*a + b'*b;
    end
    
    r = 0.5*alpha*r + 3*N;
    
end

