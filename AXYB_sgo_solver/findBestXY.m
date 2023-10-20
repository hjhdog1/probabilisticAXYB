function [X,Y,J] = findBestXY(XYset)

    N = length(XYset);
    
    Idx_min = 1;
    J_min = XYset(1).J;
    
    for i = 1:N
        if XYset(i).J < J_min
            Idx_min = i;
        end
    end
    
    X = XYset(Idx_min).X;
    Y = XYset(Idx_min).Y;
    J = XYset(Idx_min).J;

end

