function cry__WritePOSCAR(filename, crystal, varargin)
% write atomic structure defined in 'crystal' data structure to a VASP
% POSCAR file; the atomic species are sorted in ascending atomic number,
% MIND THIS when providing the POTCAR file
%
% USAGE: cry__WritePOSCAR(filename, crystal, [SDMatrix])
%
% atomic structure defined by
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.atomnum(atomID)       = atomic number of the atom 
%
% SDmatrix(atomID, 1..3) : matrix defining (and switching on) selective dynamics
%                        = 0 : coordinate is FIXED (F)
%                        = 1 : coordinate is allowed to change (T)
% 
% the atomic order for the SDmatrix corresponds to the order in
% crystal.atompos and crystal.atomnum, i.e. it is automatically resorted
% for the POSCAR 
%
% LIMITATION: number elements currently limited to 10, can be easily changed

SDmatrix = [];
switch(nargin)    
    case(2)   % default
        
    case(3)
        SDmatrix = varargin{1};
        if size(crystal.atompos) ~= size(SDmatrix)
            error('The dimension of the SDMatrix does match the dimension of the atomic positions.')
        end

    otherwise
        error('Wrong number of function arguments.')
end

%
%[atomicnum, elementstr, numatoms, atompos_sort] = getPOSCARInfos(crystal, SDmatrix);
[atomicnum, elementstr, numatoms, atompos_sort, atomnum_sort, SDmatrix_sort] = cry__GetPOSCARInfos(crystal, SDmatrix);

display(sprintf('\n   Writing to file: %s', filename))

% open file
fh = fopen(filename, 'w');

% write header plus scaling factor
fprintf(fh, '# generated with cry__WritePOSCAR\n  1.0\n');
% write supercell lattice vectors 
fprintf(fh, '  %f %f %f\n', crystal.latt'); % matrix needs to be transposed!
fprintf(fh, ' %s\n', elementstr); 
% write number of elements: LIMITED TO 10 ELEMENTS, easy to change
fprintf(fh, ' %d %d %d %d %d %d %d %d %d %d', numatoms);
if ~isempty(SDmatrix)
    fprintf(fh, '\nSelective Dynamics');
end
fprintf(fh, '\nCartesian');

% write atomic positions 
natoms = size(atompos_sort,1);
for at = 1:natoms
    % write atomic positions
    fprintf(fh, '\n %f %f %f', atompos_sort(at,:) );   
    
    % add switches for Selective Dynamics
    if ~isempty(SDmatrix_sort)
        for comp = 1:3
            if SDmatrix_sort(at, comp) == 0
                fprintf(fh, ' F ');      
            else
                fprintf(fh, ' T ');      
            end
        end
    end
end;


fclose(fh);




