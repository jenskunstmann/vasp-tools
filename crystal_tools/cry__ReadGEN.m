function crystal = cry__ReadGEN(filename)
% read in DFTB+ gen file and return its structure in the 'crystal' data
% structure
%
% filename = of gen file


% at_str2num(atom type string, atomic number) = cell array mapping the
%                               atomic symbol to the atomic number, could
%                               in principle be a generic array containing
%                               all atoms
%                            

% use Felix' routine to get the data
out = gen2mat(filename);

% get array dimensions
[natoms, three] = size(out.AtomCoordinates);

% find out which type each atom is and generate the atomnum array
atomnum = zeros(natoms,1);
for at = 1:natoms  
        atomnum(at) = cry__getAtomicNumber(out.AtomTypes(at));        
end

% save things
crystal.origin  = out.origin;
crystal.latt    = out.LatticeVectors;
crystal.atompos = out.AtomCoordinates;
crystal.atomnum = atomnum;



