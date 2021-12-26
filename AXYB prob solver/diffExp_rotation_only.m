function dexpu = diffExp_rotation_only(u)

norm_u = norm(u);
if norm_u < 1e-14
    dexpu = eye(3);
else
    norm_u_half = 0.5*norm_u;
    one_over_norm_u_sqr = 1/(norm_u^2);
    
    s = sin(norm_u_half)/norm_u_half;
    c = cos(norm_u_half);
    
    alpha = s*c;
    beta = s^2;
    
%     skew_u = so3(u);
%     dexpu = eye(3) + 0.5*beta*skew_u + one_over_norm_u_sqr *(1-alpha)*skew_u^2;

%     skew_u = [0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
%     dexpu = eye(3) + (0.5*beta)*skew_u + (one_over_norm_u_sqr *(1-alpha))*skew_u*skew_u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    skew_u = [0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
%     dexpu0 = eye(3) + (0.5*beta)*skew_u + (one_over_norm_u_sqr *(1-alpha))*skew_u*skew_u;

    dexpu = (one_over_norm_u_sqr *(1-alpha))*skew_u;
        
%     temp = 0.5*beta;
%     dexpu(1,1) = dexpu(1,1) + temp;
%     dexpu(2,2) = dexpu(2,2) + temp;
%     dexpu(3,3) = dexpu(3,3) + temp;
    
%     dexpu = dexpu*skew_u;
    dexpu = dexpu*skew_u + (0.5*beta) * skew_u;
    
    dexpu(1,1) = dexpu(1,1) + 1;
    dexpu(2,2) = dexpu(2,2) + 1;
    dexpu(3,3) = dexpu(3,3) + 1;
        
end

end