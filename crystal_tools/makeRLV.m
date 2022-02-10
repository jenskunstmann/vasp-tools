function rlv = makeRLV(lv)
% make the reciprocal basis lattice vectors (RLV)

% the vectors are given in units of the lattice parameter A
a1 = lv(1,:);
a2 = lv(2,:);
a3 = lv(3,:);

% the rlv are given in units of 2*Pi/A
% therefore we have no factor of 2*Pi
denom = dot(a1, cross(a2,a3));
rlv(1,:) = cross(a2,a3)/denom;
rlv(2,:) = cross(a3,a1)/denom;
rlv(3,:) = cross(a1,a2)/denom;

