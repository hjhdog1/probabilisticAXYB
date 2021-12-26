function [X0,Y0] = findInitGuessUsingAXXB2(A,B)
    N = size(A,3);
    U = zeros(3,3);
    V = zeros(3,3);
        
    for i = 1:N-1
        a = so3(A(:,:,1)'*A(:,:,i+1));
        b = so3(B(:,:,1)'*B(:,:,i+1));
        U = U + b*a';
        
        a = so3(A(:,:,1)*A(:,:,i+1)');
        b = so3(B(:,:,1)*B(:,:,i+1)');
        V = V + b*a';
    end
    
    X0 = eye(3);
    if det(U) > 0
        X0 = (U'*U)^(-0.5) * U';
    else
        [A, ~, B] = svd(U);        
        X0 = B * diag([1,1,-1]) * A';
    end
        
    Y0 = eye(3);
    if det(V) > 0
        Y0 = (V'*V)^(-0.5) * V';
    else
        [A, ~, B] = svd(V);        
        Y0 = B * diag([1,1,-1]) * A';
    end

end

