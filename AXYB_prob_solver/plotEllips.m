function h = plotEllips(m, axes, axes_length, lineWidth)

    if nargin < 4
        lineWidth = 1;
    end
    

    t = linspace(0, 2*pi, 100);
    cs = [cos(t); sin(t)];


    p = axes * diag(axes_length) * cs + m;
    h = plot(p(1,:), p(2,:), 'LineWidth', lineWidth);
    
%     plot(m(1), m(2), 'o');

  
  
end
