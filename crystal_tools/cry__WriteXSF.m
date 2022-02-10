function cry__WriteXSF(filename, crystal)
% write atomic structure defined in 'crystal' data structure to a xsf file
%
% atomic structure defined by
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.atomnum(atomID)       = atomic number of the atom 

display(sprintf('\n   Writing to file: %s', filename))

% open file
fh = fopen(filename, 'w');

% write header
fprintf(fh, 'CRYSTAL\n');

% write supercell lattice vectors as primitive lattice vectors
fprintf(fh, 'PRIMVEC\n');
fprintf(fh, '  %f %f %f\n', crystal.latt');

natoms = size(crystal.atompos,1);
% write atomic positions + phonon vectors
fprintf(fh, 'PRIMCOORD\n %d 1\n', natoms);
for at = 1:natoms
    fprintf(fh, '%d %f %f %f\n', crystal.atomnum(at), ...
                crystal.atompos(at,:) );           
end;


fclose(fh);
