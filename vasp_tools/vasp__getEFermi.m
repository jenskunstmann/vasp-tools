function efermi = vasp__getEFermi(filename)
% read the DOSCAR file of VASP and return the Fermi energy

% open file
fh = fopen(filename, 'r');

% skip 5 lines of header
skipline(fh, 5);

% read dimension line 
[Q,l] = fscanf(fh, ' %g  %g  %i %g %g ', [5 1]);
if l ~= 5
    error('error reading file');         
end;
%emax  = Q(1)   % maximal energy
%emin  = Q(2)   % minimal energy
%npnt  = Q(3)   % number of energy points
efermi = Q(4);  % Fermi energy

% close file
fclose(fh);
