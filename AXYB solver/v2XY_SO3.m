function [X,Y] = v2XY_SO3(v)

    X_t = zeros(3);
    Y_t = zeros(3);

    for i = 1:9
        X_t(i) = v(i);
        Y_t(i) = v(i+9);
    end
    
    X = X_t';
    Y = Y_t';

end

