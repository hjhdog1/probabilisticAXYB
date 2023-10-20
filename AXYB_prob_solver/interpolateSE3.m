function X = interpolateSE3(X1, X2, t)

% interpolate between X1 and X2 s.t. X(t=0) = X1 and X(t=1) = X2

X = eye(4);

w = LogSO3(X1(1:3,1:3)'*X2(1:3,1:3));
X(1:3,1:3) = X1(1:3,1:3) * LargeSO3(w*t);
X(1:3,4) = (1-t)*X1(1:3,4) + t*X2(1:3,4);

end

