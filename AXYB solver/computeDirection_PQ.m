function [d g H] = computeDirection_PQ(lambda,P,Q,P0,Q0,X,Y,isNewtonMethod)

    a = zeros(3,18);
    b = zeros(3,18);
    M = zeros(3,3,18);
    N = zeros(3,3,18);
    
    for i = 1:18
        a(:,i) = skew2so3(X'*P(:,:,i)' - P(:,:,i)*X);
        b(:,i) = skew2so3(Y'*Q(:,:,i)' - Q(:,:,i)*Y);
        M(:,:,i) = 0.5*(X'*P(:,:,i)' + P(:,:,i)*X);
        N(:,:,i) = 0.5*(Y'*Q(:,:,i)' + Q(:,:,i)*Y);
    end
    a0 = skew2so3(X'*P0' - P0*X);
    b0 = skew2so3(Y'*Q0' - Q0*Y);
    M0 = 0.5*(X'*P0' + P0*X);
    N0 = 0.5*(Y'*Q0' + Q0*Y);
    
    % gradient
    g = zeros(6,1);
    for i = 1:18
        g = g + lambda(i) * trace(M(:,:,i) + N(:,:,i)) * [a(:,i);b(:,i)];
    end
    g = g + [a0;b0];
    d = -g;
    
    % hessian    
    H = zeros(6);
    if isNewtonMethod        
        for i = 1:18
            H(1:3,1:3) = H(1:3,1:3) + lambda(i) * ( a(:,i)*a(:,i)' + trace(M(:,:,i) + N(:,:,i))*(M(:,:,i) - trace(M(:,:,i))*eye(3)) );
            H(4:6,4:6) = H(4:6,4:6) + lambda(i) * ( b(:,i)*b(:,i)' + trace(M(:,:,i) + N(:,:,i))*(N(:,:,i) - trace(N(:,:,i))*eye(3)) );
            H(1:3,4:6) = H(1:3,4:6) + lambda(i) * a(:,i)*b(:,i)';
        end
        H(1:3,1:3) = H(1:3,1:3) + M0 - trace(M0)*eye(3);
        H(4:6,4:6) = H(4:6,4:6) + N0 - trace(N0)*eye(3);
        H(4:6,1:3) = H(1:3,4:6)';
        
        d = -pinv(H)*g;        
    end

end

