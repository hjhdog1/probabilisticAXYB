function invT = invSE3( T )

invT = eye(4);
invT(1:3,1:3) = T(1:3,1:3)';
invT(1:3,4) = -invT(1:3,1:3)*T(1:3,4);


end

