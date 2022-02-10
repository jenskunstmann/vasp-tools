function cry__WriteGEN(filename, crystal)
% write atomic structure defined in 'crystal' data structure to a DFTB+ gen file
%
% atomic structure defined by
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.atomnum(atomID)       = atomic number of the atom 

display(sprintf('\n   Writing to file: %s', filename))

%%% get required information
natoms = size(crystal.atompos,1);

% atomicnum(1:nelements) = atomic number of the elements in the system
[ElementToAtomnum, ~ , AtomIDToElement] = unique(crystal.atomnum);
nelements = length(ElementToAtomnum);

% elementstr = string with the atomic symbols of the 'nelements' elements
elementstr = '';
for elem =  1:nelements
    elementstr = [elementstr ' ' cry__getAtomicSymbol(ElementToAtomnum(elem))];
end

% for the lines with the Cartesian coordinates
coordmat = [
    reshape(1:natoms, natoms, 1), reshape(AtomIDToElement, natoms, 1), crystal.atompos];


% problem in the code below is that every new line introduces a new line in
% the matrix -> write all in one line instead!
% if natoms == 1
%     coordmat = [(1:natoms)
%                 AtomIDToElement'    % prime unset for Nat > 1
%                 crystal.atompos' ];    
% else
%     coordmat = [(1:natoms)
%                 AtomIDToElement    % prime unset for Nat > 1
%                 crystal.atompos' ]
% end

% open file
fh = fopen(filename, 'w');

% write it out
fprintf(fh, '%d S\n', natoms);          % supercell coordinates in Cartesian coordinates
fprintf(fh, '%s\n', elementstr);        % element symbols
fprintf(fh, '%d %d %f %f %f\n', coordmat');       % Cartesian coordinates
fprintf(fh, '%f %f %f\n', crystal.origin);    % unit cell origin
fprintf(fh, '%f %f %f\n', crystal.latt');     % unit cell vectors

% close file
fclose(fh);
