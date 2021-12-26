function J = computeCost_SE3_PQ(lambda,P,Q,P0,Q0,c,X,Y)
% X,Y : SO(3)

    J = 0;
    for i = 1:18
        J = J + 0.5 * lambda(i) * ( trace(P(:,:,i)*X) + trace(Q(:,:,i)*Y) )^2;
    end
    J = J + trace(P0*X) + trace(Q0*Y) + c;

end

