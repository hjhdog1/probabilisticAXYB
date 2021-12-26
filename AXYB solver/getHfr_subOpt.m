function [H_subOpt,f_subOpt,r_subOpt] = getHfr_subOpt(H,f,r)

    % build H_subOpt
    H11 = H(1:18,1:18);
    H12 = H(1:18,19:24);
    H22 = H(19:24,19:24);
    invH22 = pinv(H22);
    
    H_subOpt = H11 - H12*invH22*H12';
    
    % build f_subOpt
    f1 = f(1:18);
    f2 = f(19:24);
    
    f_subOpt = f1 - H12*invH22*f2;
    
        
    % build r_subOpt
    r_subOpt = r - 0.5*f2'*invH22*f2;
    

end

