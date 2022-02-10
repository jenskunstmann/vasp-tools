function nntable = cry__FindNearestNeighbors(cell, cutoff, varargin)
% find all the nearest neigbors of the atoms within a sphere around each
% atom; a supercell is generated to find all the NN, its size is determined
% automatically from 'cutoff', in a pessimistic estimate = supercell is
% overstimated
%
% USAGE: nntable = cry__FindNearestNeigbors(cell, cutoff [, dim]])
%
% cell           = 'crystal' structure
% cutoff         = sphere radius to search for NN
% dim() = array defining the dimension(s) to consider, default = 1:3
%
% nntable{1..natoms}.atomnum(:)       = vector with the atomic numbers of the
%                                       neighbors
% nntable{1..natoms}.distvec(:, comp) = corresponding distance vector
%                                       pointing from the atom to the neighbors
% nntable{1..natoms}.atomID(:)        = corresponding atom ID in the
%                                       primitve cell


switch(nargin)
    case(2)
        dim = 1:3;
        
    case(3)     
        dim = varargin{3}; 
        
    otherwise
        error('Wrong number of arguments')        
end


% generate a supercell that can accommodate the cutoff radius
% first: determine the size of the supercell
latt_cell = cry__getLattConst(cell);   % get lengths of the lattice vectors

% THE FOLLOWING LINES ARE NOT SAVE TO USE
%relcoord = cell.atompos/cell.latt;    % atom position in relative coordinates
% guessed empirical correction to make supercell bigger, method below only
% works correctly for rectangular cells, for non-rectangular ones some
% neigbors were missing
%increase = 1.5;                       

% this is a very conservative estimate, assuming that atoms might sit very
% close to the cell boundaries in all directions
for i=1:3   % loop over lattice vectors    
    maxlen = ceil(cutoff/latt_cell(i));
    cmin(i) = -maxlen;
    cmax(i) = maxlen;
%     cmin(i) = floor(min(relcoord(:,i)) - cutoff*increase/latt_cell(i));
%     cmax(i) = ceil(max(relcoord(:,i)) + cutoff*increase/latt_cell(i)) - 1;
end
cmin
cmax
% now make the supercell
supercell = cry__MakeSuperCell(cell, cmin, cmax);
%cry__PlotAtomicStructure(supercell)

natoms = length(cell.atomnum);
for at = 1:natoms
    [atoms, dist_array, dist_vecs] = cry__FindAtomsNearPoint(supercell, cell.atompos(at,:), cutoff, dim);
    nnIDs = find(dist_array);    % indices of neigbors in 'atoms', remove atom itself
    nntable{at}.atomnum = supercell.atomnum(atoms(nnIDs));
    nntable{at}.distvec = dist_vecs(nnIDs,:);
    nntable{at}.atomID  = cry__FindAtomIDsInPrimitiveCell(atoms(nnIDs), natoms);   
end


% % check if supercell is big enough to accommodate the cutoff
% minlatt_cell      = min(cry__getLattConst(cell));
% minlatt_supercell = min(cry__getLattConst(supercell));
% if (2*cutoff + minlatt_cell) >= minlatt_supercell
%     error('Caution: The supercell is not big enough to find all nearest neigbors. Make it bigger or decrease the cutoff radius.')
% end


