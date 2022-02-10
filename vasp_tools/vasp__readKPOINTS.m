function [nkpnts] = vasp__readKPOINTS(filename)
% read VASP KPOINTS file
%
% return the number of k-points
%
% LIMITATIONS: currently it can only read the number of k-points from the
% second line

% open file
fh = fopen(filename, 'r');

% skip 1 lines of header
skipline(fh, 1);

% read line
[nkpnts,l] = fscanf(fh, ' %d ', 1);
if l ~= 1
    error('error reading file: number of k-points');         
end;

% close file
fclose(fh);
