function dexpu_inv = diffExpInv_rotation_only(u)

norm_u = norm(u);
if norm_u < 1e-14
    dexpu_inv = eye(3);
else
    norm_u_half = 0.5*norm_u;
    one_over_norm_u_sqr = 1/(norm_u^2);
    
    s = sin(norm_u_half)/norm_u_half;
    c = cos(norm_u_half);
    
    gamma = c/s;
    
    
    skew_u = [0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
    
    dexpu_inv = (one_over_norm_u_sqr * (1-gamma)) * skew_u;
    dexpu_inv = dexpu_inv*skew_u - 0.5*skew_u;
    
    dexpu_inv(1,1) = dexpu_inv(1,1) + 1;
    dexpu_inv(2,2) = dexpu_inv(2,2) + 1;
    dexpu_inv(3,3) = dexpu_inv(3,3) + 1;
        
end

end