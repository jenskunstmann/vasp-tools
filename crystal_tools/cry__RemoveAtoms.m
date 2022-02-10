function newpatch = cry__RemoveAtoms(oldpatch, atomIDs)
% remove atoms specified by their atomID
% atomID(:) = vector with atomic numbers
%
% atomic structure defined by
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.atomnum(atomID)       = atomic number of the atom 

newpatch = oldpatch;
newpatch.atompos(atomIDs,:) = [];
newpatch.atomnum(atomIDs) = [];
end
