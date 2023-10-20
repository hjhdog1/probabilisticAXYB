function [dexpu, C_halfv] = diffExp(u,v)

norm_u = norm(u);
% skew_v = so3(v);
skew_v = [0,-v(3),v(2); v(3),0,-v(1); -v(2),v(1),0];
if norm_u < 1e-14
    dexpu = eye(3);
    C_halfv = 0.5*skew_v;
else
    norm_u_half = 0.5*norm_u;
    one_over_norm_u_sqr = 1/(norm_u^2);
    
    s = sin(norm_u_half)/norm_u_half;
    c = cos(norm_u_half);
    
    alpha = s*c;
    beta = s^2;
    
    skew_u = so3(u);
    dexpu = eye(3) + 0.5*beta*skew_u + one_over_norm_u_sqr *(1-alpha)*skew_u^2;
    
    C_halfv = 0.5*(beta)*skew_v + (1-alpha)*one_over_norm_u_sqr * (skew_v*skew_u + skew_u*skew_v)...
              + (alpha-beta)*one_over_norm_u_sqr*(u'*v)*skew_u ...
              + one_over_norm_u_sqr*(0.5*beta - 3*(1-alpha)*one_over_norm_u_sqr)*(u'*v)*skew_u^2;
    
end

end

