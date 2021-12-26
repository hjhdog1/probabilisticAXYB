function [A,B] = randomSorting(A,B)

N = size(A,3);
x = rand(N,1);
[x, order] = sort(x);

A = A(:,:,order);
B = B(:,:,order);


end

