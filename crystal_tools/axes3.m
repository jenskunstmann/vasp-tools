function axes3()
% draws arrows that indicate the origin of the Cartesian axes

% specify the variables
ar_len = 1.5;
ar_bodywidth = 0.1;
ar_headlen = ar_len*5/12;
ar_headwidth = 2.5*ar_bodywidth;
color = 'black';

% plot the actual arrows
ar = arrow3([0 0 0], [ar_len 0 0], ar_bodywidth, ar_headlen, ar_headwidth);
set(ar, 'FaceColor', color);
ar = arrow3([0 0 0], [0 ar_len 0], ar_bodywidth, ar_headlen, ar_headwidth);
set(ar, 'FaceColor', color);
ar = arrow3([0 0 0], [0 0 ar_len], ar_bodywidth, ar_headlen, ar_headwidth);
set(ar, 'FaceColor', color);

% plot the x,y,z axes labels
shift = -ar_bodywidth;
coord = [[ar_len shift shift]' ...
         [shift ar_len shift]' ...
         [shift shift ar_len]'];
labstr = {'x' 'y' 'z'};

% the 3 rows of 'coord' are the 3 components
text(coord(1,:), coord(2,:), coord(3,:), labstr, ...
    'FontSize',20, 'FontName','Arial');