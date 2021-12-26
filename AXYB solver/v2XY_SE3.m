function [X,Y] = v2XY_SE3(v)

    X_t = zeros(3);
    Y_t = zeros(3);

    for i = 1:9
        X_t(i) = v(i);
        Y_t(i) = v(i+9);
    end
    
    X = [X_t' v(19:21,1); zeros(1,3) 1];
    Y = [Y_t' v(22:24,1); zeros(1,3) 1];

end

