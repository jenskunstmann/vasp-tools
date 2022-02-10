function out = gen2mat(file)
% read in a DFTB geoometry file
%
% out.AtomTypes{atomID} = cell array containing a string with the atomic symbol
% out.AtomCoordinates(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% out.origin(:)            = Cartesian coordinates for origin of the cell
% out.LatticeVectors(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; 
%                                    comp=1,2,3 = Cartesian component
% out.Info                 = string containing infos
% out.Name                 = filename string, that was read
%
% written by Felix Zoegiebel

fid=fopen(file);

% read number of atoms
N=fscanf(fid,'%d',1);

% read file type (F or S or C)
FSC=fgetl(fid);
FSC=FSC(FSC~=' ');

switch FSC
    case 'F'
        relcoord=true;  lattice=true;
    case 'S'
        relcoord=false; lattice=true;
    case 'C'
        relcoord=false; lattice=false;
    otherwise
        error('File-Type (character in first line) must be F, S or C. Found %s\n',FSC);
end

% read next non-empty line, which holds the atom types
s=[]; while isempty(s), s=fgetl(fid); end

% "parse" atom types
space=' ';
types={};
a=find(s~=space,1);
while ~isempty(a)
    s=s(a:end);
    b=find((s==space),1);
    if isempty(b), b=length(s); end
    types=[types {s(1:b)}];
    s=s(b+1:end);
    a=find(s~=space,1);
end

% read atom type numbers and atom positions
tXYZ = fscanf(fid,'%*d %d %g %g %g',[4 N])';

% translate atom type numbers to atom types from 2nd line
% JENS K.: strtrim introduced, DATE 5.6.13
t = strtrim(types(tXYZ(:,1)));

% extract positions
XYZ = tXYZ(:,2:4);

RepVecs=zeros(3);
if lattice
    % read origin of unit cell
    CoUC=fscanf(fid,'%g',3)';
    % read repetition vectors in rows (as in file)
    RepVecs=fscanf(fid,'%g',[3 3])';
    if relcoord
        XYZ=(RepVecs'*XYZ')';
    end
    XYZ=XYZ; % -CoUC(ones(1,N),:);
end

fclose(fid);

out.AtomTypes=t;
out.AtomCoordinates=XYZ;
out.LatticeVectors=RepVecs;
out.Info='All units in Angstrom, vectors in rows';
out.Name=file(1:end-4);
% added by Jens Kunstmann, 10.04.13
if exist('CoUC', 'var')
    out.origin=CoUC;
else
    out.origin=[0 0 0];
end

end