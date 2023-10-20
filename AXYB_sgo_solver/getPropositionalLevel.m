function level = getPropositionalLevel(Val,proposition)
% Val should be increasingly sorted.

    n = length(Val);
    if n == 1
        level = Val(1);
        return;
    end
    
    if proposition <= 0
        level = Val(1);
        return;
    end
    
    if proposition >= 1
        level = Val(n);
        return;
    end
        
    ref = (n-1)*proposition+1;
    low = floor(ref);
    high = low+1;
    frac = ref - low;
    
    level = Val(low) + (Val(high) - Val(low))*frac;


end

