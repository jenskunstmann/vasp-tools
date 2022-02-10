function [const] = cry__getLattConst(crystal)
% determine lattice constants from the lattice vectors


const(1) = norm(crystal.latt(1,:),2);
const(2) = norm(crystal.latt(2,:),2);
const(3) = norm(crystal.latt(3,:),2);

end