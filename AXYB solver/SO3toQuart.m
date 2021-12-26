function q = SO3toQuart(R)
% R : SO(3)
% q : quarternion
    q = zeros(4,1);
    q(1) = 0.5*sqrt(R(1,1) + R(2,2) + R(3,3) + 1);
    
    if(abs(q(1)) > 1e-8)
        mult = 1/(4*q(1));
        q(2) = (R(3,2) - R(2,3)) * mult;
        q(3) = (R(1,3) - R(3,1)) * mult;
        q(4) = (R(2,1) - R(1,2)) * mult;
    else
        mult = 1/sqrt(R(1,2)^2 * R(1,3)^2 + R(1,2)^2 * R(2,3)^2 + R(1,3)^2 * R(2,3)^2);
        q(2) = R(1,3)*R(1,2)*mult;
        q(3) = R(1,2)*R(2,3)*mult;
        q(4) = R(1,3)*R(2,3)*mult;
    end

end

