function J = computeCost_SE3_reduced(H,f,r,X,Y)
% X,Y : SE(3)

v = XY2v_SE3(X,Y);
J = 0.5*v'*H*v + f'*v + r;

end

