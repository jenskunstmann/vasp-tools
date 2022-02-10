function xsfconvert
% convert fractional (reduced) coordinates into Cartesian coordinates
%
% one can convert a cif file by using the Bilbao crystallographic server
% http://www.cryst.ehu.es/cryst/wpassign.html
% this tool returns fractional coordinates
% this programs converts them into Cartesian coordinates

% lattice constants
A = 4.982;
C = 10.0;

% hexagonal lattice vectors
latt(1,:) = A*[sqrt(3)/2, 0.5, 0];
latt(2,:) = A*[-sqrt(3)/2, 0.5, 0];
latt(3,:) = C*[0, 0, 1];

latt

myangle(latt(1,:),latt(2,:)) % angle between a1 and a2

% fractional coordinates
fracpos = ...
[ 0.33333,0.66667,0.500000; ...
  0.66667,0.33333,0.500000;...
  0.33066,0.00000,0.500000;...
  0.00000,0.33066,0.500000;...
  0.66934,0.66934,0.500000;...
  0.66934,0.00000,0.500000;...
  0.00000,0.66934,0.500000;...
  0.33066,0.33066,0.500000];

natoms = size(fracpos,1);   % number of atoms
cartpos = zeros(natoms,3);  % initialization

% generate Cartesian coordinates  
for at = 1:natoms
    cartpos(at,:) = fracpos(at,1)*latt(1,:) + fracpos(at,2)*latt(2,:) + ...
                    fracpos(at,3)*latt(3,:);
end;


cartpos
