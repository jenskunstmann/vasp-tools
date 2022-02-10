function newpatch = cry__AddAtom(oldpatch, atompos, atomicnum)
% add an atom to 'oldpatch' at position 'position'

newpatch = oldpatch;
natoms = length(newpatch.atomnum);
newpatch.atompos(natoms+1,:) = atompos;
newpatch.atomnum(natoms+1) = atomicnum;
end
