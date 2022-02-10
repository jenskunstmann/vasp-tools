function outstruc = cry__ShiftAtomicStructure(instruc, originshift)
% this functions does
% 1. shifts the origin of the 'instruc' by 'originshift'
% 2. sets the origin to [0 0 0] and shifts all atoms accordingly
% 3. translates all atoms back into the cell

% initialize the output and modify it 
outstruc = instruc;
outstruc.origin = [0 0 0];

natoms = size(instruc.atompos,1);
% do the shifts of steps 1 and 2
cartpos = instruc.atompos - originshift(ones(natoms,1),:) - instruc.origin(ones(natoms,1),:);
% convert in to relative coordinates
relcoord = cartpos/outstruc.latt;
% translate atoms outside of the cell back into it
relcoord = mod(relcoord,1);
% convert back into Cartesian coordinates
outstruc.atompos = relcoord * outstruc.latt;    

end
