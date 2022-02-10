function atomIDs = cry__FindAtomsByType(crystal, atomicnum)
% find all atoms of atomic number 'atomicnum'

natoms = length(crystal.atomnum);

atomIDs = [];
atomIDs_count = 0;
for nat = 1:natoms   
    if crystal.atomnum(nat) == atomicnum    
        atomIDs_count = atomIDs_count + 1;
        atomIDs(atomIDs_count) = nat;
    end
end