function invA = invertData(A)

    N = size(A,3);
    
    invA = zeros(4,4,N);
    for i = 1:N
        invA(:,:,i) = invSE3(A(:,:,i));
    end

end

