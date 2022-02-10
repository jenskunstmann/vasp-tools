function vasp_PlotProjBands_spallatomsintegrated()
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
  
% plot only a certain part of the band plot
% define starting and end points
% pnt_min = 3;
% pnt_max = npnts/2;
% % select the corresponding data points
% bands.klabels = {'\Gamma','X'};                 % special point labels
% bands.kpnt_pos = kpnt_pos([pnt_min:pnt_max],:); % positions of the kpoints   
% bands.eval = eval(pnt_min:pnt_max,:)-efermi;    % eigen energies with Fermi level at E=0
% bandchar = bandchar(pnt_min:pnt_max,:,:,:);
     
% define the atoms sum over
satoms = 1:natoms;     
     

% sigma states (s+px+py) of a single atom
bands.bchar = 0;
for atom = satoms
    for orbital=[2 3 4]
        bands.bchar = bands.bchar + bandchar(:,:,atom,orbital);
    end
end
bands.charcol = [1 0 0];
bands.charscal = .5;               % fatness of bands  = charscal * band character  
vasp__plotBands(bands);

hold on

% plot the pi bands
bands.bchar = 0;
for atom = satoms
    for orbital=[5 6 7 8 9]
        bands.bchar = bands.bchar + bandchar(:,:,atom,orbital);
    end    
end
bands.charcol = [0 0 1]; %[1 0.6 0];
bands.charscal = .2;               % fatness of bands  = charscal * band character  
virtax = vasp__plotBands(bands);  % plot bands and return the virtual axis

% plot band numbers
%vasp__PlotBandNumbers(bands, virtax)

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
title_str = sprintf('%s, \\sigma and \\pi states',SYS.ID);
title(title_str); %,'Interpreter','none');

