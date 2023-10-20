function [XSE3,YSE3] = getOptSE3(XSO3,YSO3,H,f)

    xy = -pinv(H(19:24,19:24)) * (H(19:24,1:18)*XY2v_SO3(XSO3,YSO3) + f(19:24));
    XSE3 = [XSO3 xy(1:3);zeros(1,3) 1];
    YSE3 = [YSO3 xy(4:6);zeros(1,3) 1];
    
end

