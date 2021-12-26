function v = XY2v_SE3(X,Y)

    v = zeros(24,1);
    
    X_t = X(1:3,1:3)';
    Y_t = Y(1:3,1:3)';
    
    for i = 1:9
        v(i) = X_t(i);
        v(i+9) = Y_t(i);
    end
    v(19:21,1) = X(1:3,4);
    v(22:24,1) = Y(1:3,4);
    
end

