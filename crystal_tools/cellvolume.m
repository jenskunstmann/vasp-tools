function volume = cellvolume(latt)
% calculate the volume of the unit cell given by 'latt'

volume = dot(cross(latt(1,:), latt(2,:)), latt(3,:));
