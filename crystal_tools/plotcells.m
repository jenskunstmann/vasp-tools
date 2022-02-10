function plotcells(latt, origin, cmin, cmax, color)
% plot a the unit cells vectors, defined in the crystal structure, and
% multuiples of it
%
% latt(num, comp) = lattice vector matrix, num = 1,2,3, comp - 1,2,3 ...
% Cartesian component
% origin = [x y z], vector of the origin
% cmin = [i1_start i2_start i3_start], supercells, range of plotting
% cmax = [i1_end   i2_end   i3_end]

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


hold on
for i1 = cmin(1):cmax(1)     % loop over the lattice vectors
    for i2 = cmin(2):cmax(2)
        for i3 = cmin(3):cmax(3)
            R   = origin + [i1 i2 i3]*latt;   
            
            xpos = R(1) + pnts(1,:);
            ypos = R(2) + pnts(2,:);
            zpos = R(3) + pnts(3,:);     

            % conneted the vertices one by one with lines
            % the 3 rows of 'pnts' are the 3 components
            line(xpos, ypos, zpos, 'Color', color, 'LineWidth', 3);
            
        end;
    end;
end;




