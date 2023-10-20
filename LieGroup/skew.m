function output = skew(input)

assert(isequal(numel(input), 3), 'input is not a 3d vector.');

output = [    0    -input(3)  input(2);
         input(3)  0    -input(1);
        -input(2)  input(1)  0   ];

end