function cry__AnalyzeStructure(filename, type)
% plot information about the lattice, crystal and others of crystal
% structures
%
% filename = file to read in
% type = 'vasp' for POSCAR, CONTCAR, *CAR VASP files
%      = 'xsf' for XCrysDen *.xsf files
%      = 'gen' for DFTB+ *.gen files


switch(type)
    case 'xsf' 
        crystal = cry__readXSF(filename);
    case 'vasp'
        crystal = cry__readCONTCAR(filename);
    case 'gen'        
        crystal = cry__ReadGEN(filename);
end


%%%%% print out details
natoms = length(crystal.atomnum)

% dimensions
crystal.latt
[A,B,C] = getLattConst(crystal.latt);
% angle definitions follow chrystallographic conventions
alpha = vecangle(crystal.latt(2,:),crystal.latt(3,:));
beta  = vecangle(crystal.latt(1,:),crystal.latt(3,:));
gamma = vecangle(crystal.latt(1,:),crystal.latt(2,:));

display('A, B, C, alpha, beta, gamma:')
[A,B,C, alpha, beta, gamma]

% areas 
A12 = norm(cross(crystal.latt(1,:),crystal.latt(2,:)));
A13 = norm(cross(crystal.latt(1,:),crystal.latt(3,:)));
A23 = norm(cross(crystal.latt(2,:),crystal.latt(3,:)));

display('areas: A12, A13, A23:')
[A12, A13, A23]

% the rlv are given in units of 2*Pi/Ang
rlv = cry__GetRLV(crystal)
% reciprocal areas 
rA12 = norm(cross(rlv(1,:), rlv(2,:)));
rA13 = norm(cross(rlv(1,:), rlv(3,:)));
rA23 = norm(cross(rlv(2,:), rlv(3,:)));

display('reciprocal areas: A12, A13, A23:')
[rA12, rA13, rA23]

% kpoints = 2*2
% kAreaDensity = rA12*4*pi^2/kpoints % [1/Ang^2]
% Npoints = sqrt(rA12*4*pi^2/0.1) % 0.1 is the kAreaDensity for MoS2

% %Aatomic_t = 5.045139033263441/2; % atomic area of flat triangular sheet in DFT/PBE
% Aatomic_t = 7.353411164911921/3  % atomic area of flat triangular sheet determined from hexagonal sheet in DFT/PBE
% Aatomic = A12/natoms
% eta = 2*(Aatomic - Aatomic_t)/Aatomic_t/3

% bcc_lattconst = 2*crystal.latt(1,2)
% fcc_lattconst = 2*crystal.latt(1,2)

crystal.atompos


% playing with rec. lattice
kpnt = [2/3 1/3 0]*rlv

% % try to determine the radius and origin of the tubular crystal
% % define starting values for the search, i.e.
% % origin 'orig' , axis vector 'ax' and radius 'R' of the tube
% orig = [0.5 0.5 0] * crystal.latt
% orig_x_init = orig(1); 
% orig_y_init = orig(2);
% ax_phi_init = 0;    % in units of pi
% ax_theta_init = 0; 
% R_init =  4.11; 
% 
% findtube(cartpos, [orig_x_init orig_y_init ax_phi_init ax_theta_init R_init]);

% find angle between twisted bilayers
% atnum_top = [33 85];
% atnum_bottom = [15 90];
% 
% vec_top = crystal.atompos(atnum_top(2),:) - crystal.atompos(atnum_top(1),:)
% vec_bottom = crystal.atompos(atnum_bottom(2),:) - crystal.atompos(atnum_bottom(1),:)
% angle = vecangle(vec_bottom(1:2), vec_top(1:2))
% %arrow3([0 0 0], vec_top, 1, 1, 1)
% %arrow3([0 0 0], vec_bottom, 1, 1, 1)


%%%% display result
% % color/radius defined by atom type
% atomicnummax = 120; % there are actually only 118 elements
% atomicradius = 0.5*ones(atomicnummax,1); % radius for each atom type
% atomiccolor  = zeros(atomicnummax,3);    % color for each atom type
% atomicradius(1)  = .3; atomiccolor(1,:)  = [0 0 0]; % H
% atomicradius(16) = .4; atomiccolor(16,:) = [1 1 0]; % S
% atomicradius(42) = .6; atomiccolor(42,:) = [0 0 1]; % Mo
%cry__PlotAtomicStructure(crystal); %, [], [], atomicradius, atomiccolor)

end

function [A,B,C] = getLattConst(latt)
% determine lattice constants from the lattice vectors
A = norm(latt(1,:),2);
B = norm(latt(2,:),2);
C = norm(latt(3,:),2);

end