<<<<<<< HEAD
function [ h ] = drawArrow(x, y, arrow_color)
=======
function [ h ] = drawArrow(x, y)
>>>>>>> master

h = annotation('arrow');
set(h,'parent', gca, ... 
    'position', [x(1),y(1),x(2)-x(1),y(2)-y(1)], ...
<<<<<<< HEAD
    'HeadLength', 9, 'HeadWidth', 9, 'HeadStyle', 'vback1',...
    'Linewidth', 1.5,...
    'Color', arrow_color);
=======
    'HeadLength', 12, 'HeadWidth', 7, 'HeadStyle', 'cback1',...
    'Color', [0 0.2 0.2]);
>>>>>>> master

end