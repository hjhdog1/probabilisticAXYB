function [X0, Y0] = findEigenInitXY(P,Q)

    XX = P(:,:,1);
    YY = Q(:,:,1);
    if det(P(:,:,1)) < 0
        XX = -XX;
        YY = -YY;
    end
    
    X0 = maxTraceFX(XX);
    Y0 = maxTraceFX(YY);   


end
