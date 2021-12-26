function [H,f,r] = getHfr(A,B,alpha)
% J = 0.5*v'*H*v + f'*v + r
% v constains X,Y,x,y

    % build H
    H11 = buildH11(A(1:3,1:3,:),B(1:3,1:3,:));
    H12 = buildH12(A,B,alpha);
    H22 = buildH22(A,alpha);
    
    H = [ H11 H12; H12' H22];
    
    % build f
    f = buildf(A,B,alpha);
    
    % build r
    r = buildr(A,B,alpha);
    

end

