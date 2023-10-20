function [lambda,P,Q,P0,Q0,c] = getPQ(H_subOpt,f_subOpt,r_subOpt)
% J = 0.5*sum(lambda_i*(trace(PiX) + trace(QiY))^2) + trace(P0X) +
% trace(Q0Y) + c

    % P,Q, lambda
    P = zeros(3,3,18);
    Q = zeros(3,3,18);
    
    [vec, val] = sortedEig(H_subOpt);
    for i = 1:18
        [P(:,:,i), Q(:,:,i)] = v2XY_SO3(vec(:,i));
        P(:,:,i) = P(:,:,i)';
        Q(:,:,i) = Q(:,:,i)';
    end
    lambda = diag(val);
    
    % P0, Q0, c
    [P0, Q0] = v2XY_SO3(f_subOpt);
    P0 = P0';
    Q0 = Q0';
    c = r_subOpt;


end

