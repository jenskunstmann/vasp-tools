function [crystal] = cry__readCONTCAR(input, varargin)
% read in VASP CHG/CHGCAR/ELFCAR/.. files
%
% USAGE: [crystal] = cry__readCONTCAR(filename [, typeatomnum]) OR
%        [crystal] = cry__readCONTCAR(fileID   [, typeatomnum]) 
%
% input = filename = name of ***CAR file  
% input = fileID   = file handle of ***CAR file
% typeatomnum(1..ntypes) = atomic number for each atom type
%
% atomic structure defined by
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.typeatomnum(atomID)       = atomic number of the atom 
%
% current limitations: 
% no negative scaling factor

% parse input
if ischar(input)
    % open file
    fh = fopen(input, 'r');
else
    fh = input;
end
if nargin == 2
    typeatomnum = varargin{1};
else 
    typeatomnum = [];
end

% read comment line 
comment = fgetl(fh)

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
crystal.latt = latt_tmp';

% rescale lattice vectors
if scale > 0
    crystal.latt = crystal.latt * scale;
else
    error('lattice volume scaling not yet implemented.')
end

% read atom type line(s)
% [types,ntypes] = sscanf(line, ' %d %d %d %d %d %d %d %d %d %d ', 10); %
% old style expression
line = fgetl(fh);
typesstr = regexp(line, '[0-9]+', 'match'); % cell array
ntypes = size(typesstr,2); 

if ntypes < 1           
    % if reading was not successfull we have VASP 5.2 file format
    vasp52 = true;    
    
    % read atomic symbols into cell array
    typeatomsym = regexp(line, '[A-Za-z]+', 'match');  % line with atomic symbols
    
    % read next line = atom type line
    line = fgetl(fh);
    typesstr = regexp(line, '[0-9]+', 'match'); % cell array
    ntypes = size(typesstr,2);        
    if ntypes < 1 
       error('error reading file');         
    end    

else
    % if reading was successfull we have the old VASP file format
    vasp52 = false;            
end;

% convert typestr into an numerical vector
for typ = 1:ntypes
    types(typ) = str2num(cell2mat(typesstr(typ)));
end

% total number of atoms
natoms = sum(types);

% check selective dynamics 
format = fgetl(fh);
firstcar = format(1);
if (firstcar == 'S' || firstcar == 's')
    %SelectiveDynamics = true
    display('Selective Dynamcis found.')
    % read next line as format line
    format = fgetl(fh);    
    [atompos,l] = fscanf(fh, ' %f %f %f %s %s %s ', [6 natoms]);    
    if l ~= natoms*6
        error('error reading file');         
    end;
    crystal.atompos = atompos(1:3,:)';    
    

else
    %SelectiveDynamics = false
    % read the atomic positions; 
    % reading them in as a (natoms X 3) matrix does not work for some reason
    % therefore we read it in as (3 X natoms) matrix and transpose afterwards
    [atompos,l] = fscanf(fh, ' %f  %f  %f ', [3 natoms]);
    if l ~= natoms*3
        error('error reading file');         
    end;
    crystal.atompos = atompos';
end


%%%%% convert coordinates to Cartesian coordinates

% check format of coordinates and convert if necessary
firstcar = format(1);
if (firstcar == 'C' || firstcar == 'c' || firstcar == 'K' || firstcar == 'k')
    % coordinates are given in Cartesian coordinates, 
    display('Atom position given in Cartesian coordinates');
    crystal.atompos = crystal.atompos * scale;
else    
    % coordinates are given in direct coordinates; need to be converted
    display('Atom position given in direct coordinates');
    crystal.atompos = crystal.atompos * crystal.latt;
end

% if not provided by the user get the typeatomnum() vector
if isempty(typeatomnum)
    if vasp52
        % extract the information from the atomic symbols
        for typ = 1:ntypes
            typeatomnum(typ) = cry__getAtomicNumber(cell2mat(typeatomsym(typ)));
        end
    else
        % set all to hydrogen if absolutley no information is given
        typeatomnum = ones(ntypes,1);
    end
end

% make atomnum vector
maxID = 0;
for typ = 1:ntypes
    minID = maxID + 1;
    maxID = minID + types(typ) - 1;
    crystal.atomnum(minID:maxID) = typeatomnum(typ);
end

% set origin 
crystal.origin = [0 0 0];

% close file only if filename is given
if ischar(input) 
    fclose(fh);
end


% function typeatomnum = AtsymToAtnum(typeatomsym)
% % get the atomic numbers of each type from the atomic symbols string given
% % in VASP
% % if an atomic symbol cannot be matched, the corresponding atomic number
% % will be 1
% %
% % typeatomsym(1..ntypes) = cell array of atomic symbols, one for each type
% %
% % typeatomnum(1..ntypes) = atomic number for each atom type
% 
% %adata = cry__getAtom2Num();
% 
% % get cell dimensions
% %nadata = size(adata,1);
% ntypes = size(typeatomsym,2); 
% 
% % initialize, if an atom type is not found the one will be 1 = Hydrogen
% typeatomnum = zeros(ntypes,1);
% for typ = 1: ntypes
%     typeatomnum(typ) = cry__getAtomicNumber(cell2mat(typeatomsym(typ)));
%     
% %     for ia = 1:nadata
% %         if strcmp(typeatomsym(typ), adata(ia,1))
% %             typeatomnum(typ) = cell2mat(adata(ia,2));
% %             break;
% %         end
% %     end
% % 
% %     if typeatomnum(typ) == 0
% %         error(sprintf('Atomic symbool %s could not be matched.\nCheck VASP CONTCAR file or the cry__getAtom2Num() subroutine.', cell2mat(typeatomsym(typ))))
% %     end
%     
% end

