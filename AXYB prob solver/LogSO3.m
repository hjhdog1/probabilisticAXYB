function output = LogSO3(input)
% LieGroup Matlab Library [Syungkwon Ra, dearLenin@gmail.com]

%     cosTheta = (trace(input)-1)/2;
    cosTheta = (input(1,1) + input(2,2) + input(3,3) - 1) * 0.5;
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
            output = input(1:3,k);
            output(k) = output(k)+1;
            output = output/sqrt(2*(1+input(k,k))) * theta;
        else
%             output = zeros(3,1);
%             w_hat = (input-input')*(theta/(2*sin(theta)));
%             output(1,1) = w_hat(3,2);
%             output(2,1) = w_hat(1,3);
%             output(3,1) = w_hat(2,1);
%             output = [w_hat(3,2); w_hat(1,3); w_hat(2,1)];
            output = [input(3,2)-input(2,3); input(1,3)-input(3,1); input(2,1)-input(1,2)]*(theta/(2*sin(theta)));
        end
    end    

end