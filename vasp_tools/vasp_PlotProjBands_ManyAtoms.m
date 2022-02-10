function vasp_PlotProjBands_ManyAtoms()
% plot the projected bandstructure in a matrix of subplots for many
% different atoms

global SYS

% read PROCAR file
file_procar = sprintf('%s/%s',SYS.path,SYS.procar)
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
satoms = [1:natoms]; %84-GB sep3;
%satoms = [65 73 92 93]; %84-GB sep7;
%satoms = [81 91 105 106]; %84-GB sep9;

nsatoms = length(satoms);

% dimension of the (pdiml x pdimc) subplot matrix
% this is usually a square matrix
pdimc = ceil(sqrt(nsatoms));
pdiml = pdimc;

for i = 1:nsatoms    
    % select the atom
    atom = satoms(i);         

    % sum up all or select some of the orbital characters;
    % usually we plot the total character of a particular atom
    bands.bchar = 0;
    for orbital=[1:9]
        bands.bchar = bands.bchar + bandchar(:,:,atom,orbital);
    end 

    % compile the relevant data for plotting
    bands.emax = SYS.emax;                % plotting range
    bands.emin = SYS.emin;
    bands.klabels = SYS.klabels;          % special point labels
    bands.eval = eval-efermi;         % eigen energies with Fermi level at E=0
    bands.plotbands = 1;
    bands.linespec = '-k';             % 'linespec' of the lines of the bands
    bands.kpnt_pos = kpnt_pos;        % positions of the kpoints  
    bands.charscal = 0.2;               % fatness of bands  = charscal * band character
    %bands.charcol = color(orbital,:); % default color
    bands.charcol = color(2,:);        % modified color
    
    % create the subplot
    subplot(pdiml,pdimc,i);

    % plot the bands
    vasp__plotBands(bands);

    % draw title
    title_str = sprintf('%s  - atom %d',SYS.ID, atom);
    title(title_str); %,'Interpreter','none');
end

