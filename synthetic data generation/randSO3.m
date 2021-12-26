function R = randSO3()
% generate random SO(3)

% % method 1
%     w = randn(3,1);
%     w = w/norm(w);
%     t = rand()*pi;
%     w = w*t(1);
% 
%     R = LargeSO3(w);

% method 2
    t = rand()*2*pi;
    st = sin(t);
    ct = cos(t);
    
    v = randn(3,1);
    v = v/norm(v);
    
    x = [1;0;0] - v;
    x = x/norm(x);
    
    R = (2*x*x' - eye(3))*[1 0 0; 0 ct -st; 0 st ct];
    
% % method 3
%     R = randn(3,3);
%     R(:,1) = R(:,1)/norm(R(:,1));
%     R(:,2) = R(:,2) - (R(:,1)'*R(:,2))*R(:,1);
%     R(:,2) = R(:,2)/norm(R(:,2));
%     R(:,3) = cross(R(:,1),R(:,2));
    


end

