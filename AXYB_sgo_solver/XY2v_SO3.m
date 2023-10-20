function v = XY2v_SO3(X,Y)

    v = zeros(18,1);
    
    X_t = X';
    Y_t = Y';
    
    for i = 1:9
        v(i) = X_t(i);
        v(i+9) = Y_t(i);
    end
end

