function obj = plotcell(origin, latt)
% plot a cell spanned by the vectors defined by 'latt' at 
% the point 'origin'; all are row vectors
%
% latt(numvec,comp)
% origin(comp)

% get the LV as column vectors
a1 = latt(1,:)'; 
a2 = latt(2,:)';
a3 = latt(3,:)';
O  =  [0 0 0]';

% specify the vertices of the cell
% as array of column vectors, 
% i.e. 'pnts' is a matrix of dimension 3xN
pnts = [O a1 a1+a2 a2 O ...        
       a3 a3+a1 a3+a1+a2 a3+a2 a3 ...  
       a3 a3+a1 a1 a1+a2 a1+a2+a3 a2+a3 a2];


xpos = origin(1) + pnts(1,:);
ypos = origin(2) + pnts(2,:);
zpos = origin(3) + pnts(3,:);

% conneted the vertices one by one with lines
% the 3 rows of 'pnts' are the 3 components
obj = line(xpos, ypos, zpos, 'Color', 'black', 'LineWidth', 3);




