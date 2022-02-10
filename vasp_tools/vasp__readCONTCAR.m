function [structure] = vasp__readCONTCAR(input)
% read in VASP CHG/CHGCAR/ELFCAR/.. files
%
% USAGE: [structure] = vasp__readCONTCAR(filename) OR
%        [structure] = vasp__readCONTCAR(fileID) 
%
% filename = name of CHGCAR file  
%
% structure.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% structure.types(type)       = number of atoms per type 'type'
% structure.atompos(at,comp)  = atomic coordinates as given in the file, i.e. either
%                     in Cartesian or in direct coordinates
% structure.comment      = first line of the file (comment)
% structure.format       = format line, i.e. 'Direct', 'Cartesian', etc. 
% structure.vasp52       = true/false; which file format
% structure.atsymb       = line with the atomic symbols if vasp52=true
%
% current limitations: 
% only 10 atom types can be read (easy to increase)
% 'selective dynamics' is not recognized
% the atomic positions are not rescaled/converted

% parse input
if ischar(input)
    % open file
    fh = fopen(input, 'r');
else
    fh = input;
end

% read comment line 
structure.comment = fgetl(fh);

% read the scaling factor
[scale,l] = fscanf(fh, ' %f ', 1);
if l ~= 1
    error('error reading file');         
end;

% read the lattice vectors
[latt_tmp,l] = fscanf(fh, ' %f  %f  %f  ', [3 3]);
if l ~= 9
    error('error reading file');         
end;
structure.latt = latt_tmp';

% rescale lattice vectors
if scale > 0
    structure.latt = structure.latt * scale;
else
    error('lattice volume scaling not yet implemented.')
end

% read atom type line(s)
line = fgetl(fh);
[structure.types,ntypes] = sscanf(line, ' %d %d %d %d %d %d %d %d %d %d ', 10);
if ntypes < 1           
    % if reading was not successfull we have VASP 5.2 file format
    structure.vasp52 = true;    
    structure.atsymb = line;  % line with atomic symbols
    % read next line = atom type line
    line = fgetl(fh);
    [structure.types,ntypes] = sscanf(line, ' %d %d %d %d %d %d %d %d %d %d ', 10);
    if ntypes < 1 
       error('error reading file');         
    end    
else
    % if reading was successfull we have the old VASP file format
    structure.vasp52 = false;
end;
natoms = sum(structure.types);

% read format line; 'selective dynamics' will not work, yet
structure.format = fgetl(fh);

% read the atomic positions; 
% reading them in as a (natoms X 3) matrix does not work for some reason
% therefore we read it in as (3 X natoms) matrix and transpose afterwards
[atompos,l] = fscanf(fh, ' %f  %f  %f ', [3 natoms]);
if l ~= natoms*3
    error('error reading file');         
end;
structure.atompos = atompos';


% no rescaling of atomic positions so far

% close file only if filename is given
if ischar(input) 
    fclose(fh);
end


