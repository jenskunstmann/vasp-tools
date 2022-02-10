function vasp_PlotProjBands_allatomsintegrated()
% plot the projected bandstructure of one specific atom and show sigma
% (blue) and pi (orange) states as fatbands; also plot the band numbers and
% the position of specific k-points on top of the band structure
%
% also plots band numbers and kpoint positions for the selection of bands
% and kpoints for the construction of partial charge densities

global SYS

% read PROCAR file
file_procar = sprintf('%s/%s',SYS.path,SYS.procar);
[kpnt_pos, eval, bandchar, bandphase] = vasp__readPROCAR(file_procar, SYS.lorbit);

% get Fermi energy from DOSCAR
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar);
efermi = vasp__getEFermi(file_doscar);
efermi = 0

% extract data
natoms = size(bandchar,3);  % number of atoms
npnts = size(kpnt_pos,1);   % number of kpoints
nbands = size(eval,2);    % number of bands

% definition of the orbital characters:
% syntax: bandchar(kpnt,band,atom,orbital)  
% orbital = 1 (s)   -  grey
%           2 (py)  -  red
%           3 (pz)  -  green
%           4 (px)  -  blue
%           5 (dxy) -  yellow
%           6 (dyz) -  magenta
%           7 (dz2) -  cyan
%           8 (dxz) -  orange
%           9 (dx2) -  brown
color = [.5 .5 .5; 1 0 0; 0 1 0; 0 0 1; 1 1 0; ...
         1 0 1; 0 1 1; 1 0.5 0; .4 .2  0];
names = {'s', 'p_y', 'p_z', 'p_x', 'd_{xy}', 'd_{yz}', ...
         'd_{z^2}', 'd_{xz}', 'd_{x^2-y^2}'};     

% compile the relevant data for plotting     
bands.emax = SYS.emax;                % plotting range
bands.emin = SYS.emin;
bands.linespec = '-k';             % 'linespec' of the lines of the bands
bands.klabels = SYS.klabels;          % special point labels
bands.kpnt_pos = kpnt_pos;        % positions of the kpoints
bands.eval = eval-efermi;         % eigen energies with Fermi level at E=0
bands.plotbands = 1;
  


     

% all states
% long/slow version
% bands.bchar = 0;
% for atom = 1:natoms
%     for orbital=[1:9]
%         % bandchar(kpnt,band,atom,orbital) 
%         bands.bchar = bands.bchar + conj(bandphase(:,:,atom,orbital)).*bandphase(:,:,atom,orbital);
%     end
% end

% fast version of the summation
bands.bchar = sum(sum(conj(bandphase(:,:,:,:)).*bandphase(:,:,:,:),4),3);

bands.charcol = [1 0 0];
bands.charscal = .5;               % fatness of bands  = charscal * band character  
x = vasp__plotBands(bands);

% 
kpoint = 30;
band = 9;
bands.bchar(kpoint, :)

hold on
% 

% plot band numbers
%vasp__plotBandNumbers(bands, x, [0 0 0], 0)

% plot position of specific k-points on top of the band structure;
% get the position and kpoint number from the OUTCAR of the static
% calculation, 'reciprocal coordinates' = units of the RLV 
% kpos = [0 0.2 0.4];%[0.00000000 0.14285714  0.28571429 0.42857143]; %
% nkpnts = length(kpos);
% up = SYS.emax*ones(1,nkpnts);
% down = SYS.emin*ones(1,nkpnts);
% hold on
% plot([kpos;kpos], [up;down], '-r');

% draw title
title_str = sprintf('%s, all states',SYS.ID);
title(title_str); %,'Interpreter','none');

