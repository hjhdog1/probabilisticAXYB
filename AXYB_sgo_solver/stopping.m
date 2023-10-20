function isStopping = stopping( w,N )

    Ew = w*(N-1)/(N-w-2);
    
    disp(['Number of local minima found until now: ', num2str(w)]);
    disp(['Expected number of entire local minima: ', num2str(Ew)]);
    
    isStopping = 0;
    if (Ew/w) < 1.05
        isStopping = 1;
        disp('Terminating optimization: No more local minimum is expected to exist.');
    end
    
    disp(' ');
    
end

