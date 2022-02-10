function [crystal] = cry__readXSF(filename)
% read in XCrysDen xsf files
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
% - to empty lines or comments after the XSF keyword recognized
% - no atomic symbols recognized

% read file content into string
contstr = fileread(filename);
%contstrlen = length(contstr);

% check if this is a crystal structure, actually unneccessay
%returnstr = findstring(contstr, 'CRYSTAL');

% read the lattice vectors
returnstr = findstring(contstr, 'PRIMVEC');
[latt_tmp,l] = sscanf(returnstr, ' %f  %f  %f ', [3 3]);
if l ~= 9
    error('error reading file');         
end;
crystal.latt = latt_tmp';


returnstr = findstring(contstr, 'PRIMCOORD');
% read atomnum line
[A, count, ~, nextindex] = sscanf(returnstr, ' %d %d ', [1 2]);
if count ~= 2
    error('error reading file');         
end;
natoms = A(1);
% read coordinates and atomic number
[A, count] = sscanf(returnstr(nextindex:length(returnstr)), ' %d %f %f %f ', [4 natoms]);
if count ~= 4*natoms
    error('error reading file');         
end;

crystal.atomnum = cast(A(1,:),'int8');
crystal.atompos = A(2:4,:)';
crystal.origin = [0 0 0];



function returnstr = findstring(str, searchstr)
% tries to find 'searchstr' in a bigger sting 'str', returns whatever
% remains after 'searchstr' or  prints an error message if the string is
% not found 

index = strfind(str, searchstr);
if isempty(index)
    error(sprintf('XSF keyword <%s> not found.', searchstr))
end

startindex = index(1)+length(searchstr);
endindex   = length(str);
returnstr = str(startindex:endindex);
