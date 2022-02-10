function myarrow = arrow3(pnt1, pnt2, bodyrad, headlen, headrad)

% generates an arrow pointing from pnt1 to pnt2
% - for only 2 arguments it will plot an arrow 
%   with default dimensions
% - using 5 arguments the dimensions of the arrow can
%   be defined
%
% pnt1 = starting point of arrow
% pnt2 = ending point of arrow
% bodyrad = radius of the body of the arrow
% headlen = absolute length of the head part
% headrad = radius of the head
%
% myarrow = object handle

% determine the arrow direction and its (total) length
direc = pnt2 - pnt1;    
totlen = norm(direc);

% parse input arguments
switch nargin
    case 2      % use default arrow dimensions
        bodyrad = 0.01*totlen;    % radius of body
        headlen = .15*totlen;     % length of head
        headrad = 3*bodyrad;      % radius of head
    case 5
        % do nothing
    otherwise
        error('Wrong number of arguments');
end;

% 'cylinder' generates a unit rotation body 
% according to the radial function r = r(i), i = 1..npnts;
% We generate a step function and construct the arrow
npnts = 100;     % number of points in r(i)
isp = round(headlen/totlen*npnts);  % i, where the step is

r(1:isp) = headrad/isp*((1:isp)-1);     % linear rising
r(isp+1:npnts) = bodyrad;             % constant body

% generate the coordinates of an arrow 
% pointing from [0 0 1] to [0 0 0] (upside down)
[X,Y,Z] = cylinder(r,10);

% revert its direction, scale its length, 
% and shift it from the origin to 'pnt1'
X = X + pnt1(1);
Y = Y + pnt1(2);
Z = (1-Z)*totlen + pnt1(3);

hold on
% make the 3D object
myarrow = surf(X,Y,Z, 'EdgeColor','none');%, ...
%    'FaceColor','blue', ,'AmbientStrength',.2);

% rotate the arrow such that it is parallel to 'direc'
% cmd: rotate(obj, rotaxis, rotangle, origin)
% therefore we construct the rotation axis being normal to
% the z axis (where the arrow originally points to) and 
% 'direc' (where it should point to) by using the cross product
% the rotation angle is then found by the dot product 
% between both vectors
rotaxis = cross([0 0 1], direc);
angle = acos(dot([0 0 1], direc)/totlen)/pi*180;
if angle ~= 0   % otherwise rotate aborts
    rotate(myarrow, rotaxis, angle, pnt1);
end;    

hold off


