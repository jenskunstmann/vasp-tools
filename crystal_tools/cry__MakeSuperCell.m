function supercell = cry__MakeSuperCell(crystal, cmin, cmax)
% create a super cell, origin is at cmin*latt
%
% crystal = 'crystal' data structure of primitive structure
% cmin = [i_a1, i_a2, i_a3]
% cmax = [i_a1, i_a2, i_a3]
% 
% supercell = 'crystal' data structure of super structrue

% initialize data
natoms = length(crystal.atomnum);
multiplicity_vec = cmax - cmin + ones(1,3); % number of repititions in each lattice direction
multiplicity = prod(multiplicity_vec);      % number primitive cells
supercell.atompos   = zeros(natoms*multiplicity, 3);
supercell.atomnum   = zeros(natoms*multiplicity, 1);
supercell.latt(1,:) = multiplicity_vec(1) * crystal.latt(1,:);
supercell.latt(2,:) = multiplicity_vec(2) * crystal.latt(2,:);
supercell.latt(3,:) = multiplicity_vec(3) * crystal.latt(3,:);
supercell.origin    = cmin * crystal.latt;

% loop over the supercell lattice vectors
imin = 1;
for i1 = cmin(1):cmax(1)     
    for i2 = cmin(2):cmax(2)        
        for i3 = cmin(3):cmax(3)                                          
            R = [i1 i2 i3]*crystal.latt;
            supercell.atompos(imin:(imin+natoms-1),:) = R(ones(natoms,1),:) + crystal.atompos;
            supercell.atomnum(imin:(imin+natoms-1)) = crystal.atomnum;            
            imin = imin + natoms;
        end
    end
end
