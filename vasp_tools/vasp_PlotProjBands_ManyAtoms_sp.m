function vasp_PlotProjBands_ManyAtoms_sp()
% plot the projected bandstructure for many different atoms in a matrix of
% subplots; show sigma (blue) and pi (orange) states as fatbands

global SYS

% read PROCAR file
file_procar = sprintf('%s/%s',SYS.path,SYS.procar);
[kpnt_pos, eval, bandchar] = vasp__readPROCAR(file_procar, SYS.lorbit);

% get Fermi energy from DOSCAR
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar);
efermi = vasp__getEFermi(file_doscar);

natoms = size(bandchar,3);  % number of atoms

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

% define which atoms to plot
%satoms = [7,46,13,52,19,58,25,64,31,69,35,72];   % 'pureZGNR'
%satoms = [3,23,49,52,39,29,55,58,15,36];   % 'ZGNR+HT'
%satoms = [3,23,55,61,58,39,29,49,52,36];   % 'ZGNR+HT+BC'
%satoms = [35,33,29,25,21,17,13,9,5];   % 'pureAGNR'
%satoms = [68,64,60,56,52];   % 'pureAGNR'
%satoms = [62,58,54,72,73,71,45,41,37,35]; % 'AGNR+BC' x-Richtung
%satoms = [62,58,54,72,73]; % 'AGNR+BC' x-Richtung
%satoms = [64,20,50,68,73]; % 'AGNR+BC' diagonal
%satoms = [26,3,4,27,5,30,8,31,9,34,12,35];% 'ZGNR.57' y-Richtung
%satoms = [50,26,3,4,27,5,30,8,31,9,34,12,35,13];% 'ZGNR.57+H' y-Richtung
%satoms = [1,26,3,25];% 'ZGNR.57+H' x-Richtung
%satoms = [7,41,8,46,13,47,14,52,19];% 'ZGNR.211' y-Richtung
%satoms = [1,3,7];% 'ZGNR.211' x-Richtung
satoms = [1:8];% 'alpha1/res1'
nsatoms = length(satoms);

% dimension of the (pdiml x pdimc) subplot matrix
% this is usually a square matrix
pdimc = ceil(sqrt(nsatoms));
pdiml = pdimc-1;

bands.emax = SYS.emax;                % plotting range
bands.emin = SYS.emin;
bands.linespec = '-k';             % 'linespec' of the lines of the bands
bands.klabels = SYS.klabels;          % special point labels
bands.kpnt_pos = kpnt_pos;        % positions of the kpoints
bands.eval = eval-efermi;         % eigen energies with Fermi level at E=0
bands.plotbands = 1;

for i = 1:nsatoms    
    % select the atom
    atom = satoms(i);         

    % create the subplot
    subplot(pdiml,pdimc,i);
    
    % sigma states (s+px+py) of a single atom
    bands.bchar = 0;
    for orbital= [1 2 4] %for sp
        bands.bchar = bands.bchar + bandchar(:,:,atom,orbital);
    end
    bands.charcol = [1 0.5 0];
    bands.charscal = 2;               % fatness of bands  = charscal * band character  
    vasp__plotBands(bands);

    hold on

    % plot the pi bands
    orbital = 3;
    bands.bchar = bandchar(:,:,atom,orbital);
    bands.charcol = [0 1 0];
    bands.charscal = 3;               % fatness of bands  = charscal * band character  
    vasp__plotBands(bands);    

    % draw title
    title_str = sprintf('%s  - atom %d',SYS.ID, atom);
    title(title_str); %,'Interpreter','none');
end

