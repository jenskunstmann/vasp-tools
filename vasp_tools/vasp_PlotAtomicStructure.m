function vasp_PlotAtomicStructure()
% plots the atomic structure and the unit cell from the CONTCAR file with
% their atomic numbers; cell repetitions are possible, different species
% have different colors and radii, atomic position are diaplayed in Angstrom
% code largely taken over from 'kvis'
%
% if modifications to the standard plotting parameters are required, it is
% better to use the cry__PlotAtomicStructure() directly
%
% the atomic numbers can also be determined in xcrysden, 
% but ONLY IF only one unit cell is drawn

global SYS

% supercells, range of plotting
cmin = [0 0 0];
cmax = [0 0 0];

file_contcar = sprintf('%s/%s',SYS.path,SYS.contcar)
crystal = cry__readCONTCAR(file_contcar);

%%%% display result
% % color/radius defined by atom type
% atomicnummax = 120; % there are actually only 118 elements
% atomicradius = 0.5*ones(atomicnummax,1); % radius for each atom type
% atomiccolor  = zeros(atomicnummax,3);    % color for each atom type
% atomicradius(1)  = .3; atomiccolor(1,:)  = [0 0 0]; % H
% atomicradius(16) = .4; atomiccolor(16,:) = [1 1 0]; % S
% atomicradius(42) = .6; atomiccolor(42,:) = [0 0 1]; % Mo

%cry__PlotAtomicStructure(crystal, cmin, cmax, atomicradius, atomiccolor)
cry__PlotAtomicStructure(crystal, cmin, cmax); 


