function output = so3(input)
% so(3) matrix form converting from 3*1 to 3*3
% R = so3(v)
% v : 3*1, R : 3*3
% LieGroup Matlab Library [Syungkwon Ra, dearLenin@gmail.com]

if (size(input) == [3,1])
    output = [    0    -input(3)  input(2);
             input(3)  0    -input(1);
            -input(2)  input(1)  0   ];
    return;
elseif (size(input) == [3,3])
    cosTheta = (trace(input)-1)/2;
    if cosTheta > 1
        cosTheta = 1;
    else 
        if cosTheta < -1
            cosTheta = -1;
        end
    end
    theta = acos(cosTheta);
    if abs(theta) < 10^-7
        output = zeros(3,1);
    else
        if abs(theta-pi) < 10^-7
            for k = 1:3
                if abs(1+input(k,k)) > 10^-7
                    break;
                end
            end
            output = input(:,k);
            output(k) = output(k)+1;
            output = output/sqrt(2*(1+input(k,k))) * theta;
        else
            output = zeros(3,1);
            w_hat = (input-input')/(2*sin(theta))*theta;
            output(1,1) = w_hat(3,2);
            output(2,1) = w_hat(1,3);
            output(3,1) = w_hat(2,1);
        end
    end    
    return;
else error('check the size of matrix');
end;