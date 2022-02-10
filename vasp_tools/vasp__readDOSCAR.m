function [e dos idos natoms] = vasp__readDOSCAR(filename)
% read the density of states (DOS) data from an DOSCAR file of VASP
%
% e(1..npnt) : column vector with the energies, E_F = 0
%  dos(1..npnt, 1..2) : column vector(s) of the dos 
% idos(1..npnt, 1..2) : column vector(s) of the integrated dos
% for spin-polarized calcs. dos and idos have two columns
%
% the Fermi level is shifted to 0 
 

% open file
fh = fopen(filename, 'r');

% read first line and extract the number of atoms per cell
[Q,l] = fscanf(fh, ' %i  %i  %i %i ', [4 1]);
if l ~= 4
    error('error reading file');         
end;
natoms = Q(1);     % number of atoms per cell

% skip 4 lines of header
skipline(fh, 4);

% read dimension line 
[Q,l] = fscanf(fh, ' %g  %g  %i %g %g ', [5 1]);
if l ~= 5
    error('error reading file');         
end;
emax  = Q(1);   % maximal energy
emin  = Q(2);   % minimal energy
npnt  = Q(3);   % number of energy points
efermi = Q(4)  % Fermi energy

% figure out if file contains spinpolarized DOS or not
% read first data line
line = fgetl(fh);
[Q,l] = sscanf(line, ' %g %g %g %g %g ');
if l == 3           % three columns mean non-spinploarized
    spinpol = 0
elseif l == 5
    spinpol = 1    % five columns mean spinploarized
else
    error('error reading file');
end

% return to beginning of the file and skip header
fseek(fh, 0, 'bof');
skipline(fh, 6);

% now read the data
if (spinpol)        % read data for the spinpolarized case
    
    [Q,l] = fscanf(fh, ' %g %g %g %g %g ', [5 npnt]);
    if l ~= 5*npnt
        error('error reading file');                 
    end;  
    
    % extract the columns
    e = Q(1,:)' - efermi;       % set Fermi level to 0.0
    dos(:,1) = Q(2,:)/natoms;   % normalize the DOS to states/(eV atom)
    dos(:,2) = Q(3,:)/natoms; 
    idos(:,1) = Q(4,:)/natoms;  % normalize the iDOS to to states/atom
    idos(:,2) = Q(5,:)/natoms;    
    
else    % read data for non-spinpolarized case
    
    [Q,l] = fscanf(fh, ' %g %g %g ', [3 npnt]);
    if l ~= 3*npnt
        error('error reading file');                 
    end;  
        
    % extract the columns
    e    = Q(1,:)' - efermi;    % set Fermi level to 0.0
    dos  = Q(2,:)'/natoms;      % normalize the DOS to states/(eV atom)
    idos = Q(3,:)'/natoms;      % normalize the iDOS to to states/atom    
    
end

% close file
fclose(fh);

return;


