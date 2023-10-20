function [vec, val] = sortedEig(A)
    
    % for symmetric matrix only

    d = size(A,1);
    [vec, val] = eig(A);
    val = real(val);
    
    for i = 1:d-1
        for j = 1:d-i
            
            if val(j,j) > val(j+1,j+1)
                tempVec = vec(:,j);
                tempVal = val(j,j);
                
                vec(:,j) = vec(:,j+1);
                vec(:,j+1) = tempVec;
                
                val(j,j) = val(j+1,j+1);
                val(j+1,j+1) = tempVal;
            end
            
        end
    end
    
    


end

