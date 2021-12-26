function XYset_new = addXY(XYset,X,Y,J)

    XYset_new = XYset;
    N = length(XYset_new);

    % add X,Y,J if XYset is empty
    if N == 0
        XYset_new = struct;
        XYset_new(1).X = X;
        XYset_new(1).Y = Y;
        XYset_new(1).J = J;
        
        return;
    end
    
    
    % return without adding if there already exists the same minimum
    for i = 1:N
        
        wX = so3(XYset(i).X' * X);
        wY = so3(XYset(i).Y' * Y);
        if norm([wX;wY]) < 1e-1
            
            if J < XYset(i).J
                XYset_new(i).X = X;
                XYset_new(i).Y = Y;
                XYset_new(i).J = J;
            end
            
            return;
        end
        
    end
    
    % add X,Y,J
    XYset_new(N+1).X = X;
    XYset_new(N+1).Y = Y;
    XYset_new(N+1).J = J;

end

