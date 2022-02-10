function [P, eval, occ] = vasp__readWAVEDERF(filename)
% Reading VASP's formatted WAVEDER (WAVEDERF) file that contains the
% optical transition matrix elements
% see: https://cms.mpi.univie.ac.at/wiki/index.php/WAVEDER
%
% P(f_band, i_band, comp) - transition matrix elements
%                        = hbar/m <final| p_c | initial>, 
%                          with the Cartesian momentum operator p_c
%                         unit of P = [eV Ang]
%                         f_band = 1:nbands - final band number, lines
%                         i_band = 1:nbands - initial band number, columns
%                         comp = 1,2,3 - cartesian component x, y, z  
%                         P is hermitean, i.e. P_ij = (P_ji)* 
% eval(band)- band energies in eV, bands=1:nbands
% occ(band) - band occupation, bands=1:nbands
%
% HINT: only the data of the first given k-point is read in; 
% multiple k-pints per file are not tested

% open file
fh = fopen(filename, 'r');

% read first line and extract the number of atoms per cell
[Q,l] = fscanf(fh, ' %i  %i  %i  ', [3 1]);
if l ~= 3
    error('error reading file');         
end
% some other of these 3 parameters should contain the k-point number
nkpnts = 1;
nbands = Q(3);     % number of bands

% read in all data 
% Q = Q_file' that is Q(columns in file, lines in file)
[Q,l] = fscanf(fh, ' %g ', [12 nbands*nbands*nkpnts]);
if l ~= 12*nbands*nbands*nkpnts
    error('error reading file');         
end
dat = Q'; % this improves the debugging of the code

% close file
fclose(fh);


% pick the right lines from the matrix
kpnt= 1;    % k-point to consider
start = 1 + (kpnt-1)*nbands^2;  % starting index
stop   = kpnt * nbands^2;       % final index

% extract data
eval = dat(start:start-1+nbands, 5);
occ = dat(start:start-1+nbands, 6);
% construct 'vector' of complex valued transition matrix elements 
% mat(lines, comp)
mat_vec = [dat(start:stop, 7) + 1i*dat(start:stop, 8), ...
              dat(start:stop, 9) + 1i*dat(start:stop, 10), ...
              dat(start:stop, 11)+ 1i*dat(start:stop, 12)];          
     
% convert into dipole matrix elements 
% <f|P|i> = mat_vec * |E_i - E_f|
% mat_vec is hermitean but the energy difference without abs() would turn P
% into an antihermitean matrix, i.e., P_ij = -(P_ji)* ; this is incorrect
P_vec =  mat_vec.* abs(dat(start:stop, 5) - dat(start:stop, 2));  
        
% assemble into final array
P = zeros(nbands, nbands, 3);
P(:,:,1) = reshape(P_vec(:,1), nbands, nbands).';
P(:,:,2) = reshape(P_vec(:,2), nbands, nbands).';
P(:,:,3) = reshape(P_vec(:,3), nbands, nbands).';


return




