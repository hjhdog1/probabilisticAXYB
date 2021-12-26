function R = Quart2SO3(q)
% R : SO(3)
% q : quarternion

    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);
    
    R = [1 - 2*qy^2 - 2*qz^2	2*qx*qy - 2*qz*qw	2*qx*qz + 2*qy*qw;
    2*qx*qy + 2*qz*qw	1 - 2*qx^2 - 2*qz^2	2*qy*qz - 2*qx*qw;
    2*qx*qz - 2*qy*qw	2*qy*qz + 2*qx*qw	1 - 2*qx^2 - 2*qy^2];

end

