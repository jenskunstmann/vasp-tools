function vasp_PlotProjBands()
% plot the projected bandstructure of one specific atom and one specific
% state

global SYS

% read PROCAR file
file_procar = sprintf('%s/%s',SYS.path,SYS.procar)
[kpnt_pos, eval, bandchar] = vasp__readPROCAR(file_procar, SYS.lorbit);

% get Fermi energy from DOSCAR
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar);
efermi = vasp__getEFermi(file_doscar)

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
         1 0 1; 1 0.5 0; 0 1 1; .4 .2  0];
names = {'s', 'p_y', 'p_z', 'p_x', 'd_{xy}', 'd_{yz}', ...
         'd_{z^2}', 'd_{xz}', 'd_{x^2-y^2}'};     

%%% START: USER INPUT %%%%%     
%make selections for plotting:
orbital = 2;
atom = 1;
bands.bchar = bandchar(:,:,atom,orbital);
title_str = sprintf('%s - atom %d, %s character',SYS.ID, atom, names{orbital});

% integrated characters:

%sigma states (s+px+py) of a single atom
% bands.bchar = 0;
% for orbital=[1 2 4]
%     bands.bchar = bands.bchar + bandchar(:,:,atom,orbital);
% end
% title_str = sprintf('%s  - sigma bands of atom %d',SYS.ID, atom);

%orbital character of ALL atoms
%= sum up all atomic contributions   
% bands.bchar = 0;
% for atom=1:natoms
%    bands.bchar = bands.bchar + bandchar(:,:,atom,orbital);
% end

%title_str = sprintf('%s  -  %s',SYS.ID, names{orbital});

% % total contribution of single atoms
% % = sum up all orbital characters
% bands.bchar = 0;
% for orbital=1:9
%     bands.bchar = bands.bchar + bandchar(:,:,atom,orbital);
% end
% title_str = sprintf('%s  - total character of atom %d',SYS.ID, atom);


% % total contribution of a set of atoms
% % = sum up all orbital characters
% orbitals = [1:9];
% atoms = [17 19];
% bands.bchar = 0;
% for atom = atoms
%     for orbital = orbitals
%         bands.bchar = bands.bchar + bandchar(:,:,atom,orbital);
%     end
% end
% title_str = sprintf('%s  - total character of atom %d',SYS.ID, atoms);  
% % define color by resetting 'orbital' variable
% orbital = 4;
    

%%% END: USER INPUT %%%%%

% compile the relevant data for plotting
bands.emin = SYS.emin;
bands.emax = SYS.emax;                % plotting range
bands.klabels = SYS.klabels;          % special point labels
bands.eval = eval-efermi;         % eigen energies with Fermi level at E=0
bands.linespec = '-k';             % 'linespec' of the lines of the bands
bands.kpnt_pos = kpnt_pos;        % positions of the kpoints  
bands.charscal = 10;               % fatness of bands  = charscal * band character
bands.plotbands = 1;
bands.charcol = color(orbital,:);

% plot the bands
vasp__plotBands(bands);

% draw title
title(title_str); %,'Interpreter','none');
SetFontsInFigure(20);

