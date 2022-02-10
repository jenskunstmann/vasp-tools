function rlv = getRLV(latt)
% determine the reciprocal lattice vectors (RLV)

% the vectors are given in Angstroms 
a1 = latt(1,:);
a2 = latt(2,:);
a3 = latt(3,:);

% the rlv are given in units of 2*Pi/Ang
% therefore we have no factor of 2*Pi
denom = dot(a1, cross(a2,a3));
rlv(1,:) = cross(a2,a3)/denom;
rlv(2,:) = cross(a3,a1)/denom;
rlv(3,:) = cross(a1,a2)/denom;

