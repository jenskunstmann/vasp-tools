function atomIDs_pc = cry__FindAtomIDsInPrimitiveCell(atomIDs_sc, natoms_pc)
% get the equivalent atom IDs of atoms in the primitive cell 'atomIDs_pc'
% from atomIDs of atoms in the supercell 'atomIDs_sc'
%
% natoms_pc = number of atom in the primitive cell

atomIDs_pc = mod(atomIDs_sc, natoms_pc);

% correction of the mod operation
for i=1:length(atomIDs_pc)
    if atomIDs_pc(i) == 0
        atomIDs_pc(i) = natoms_pc;
    end
end