function vasp_PlotProjBands_ManyAtoms_allcharacters()
% plot the projected bandstructure in a matrix of subplots for many
% different atoms or orbitals

global SYS

% read PROCAR file
file_procar = sprintf('%s/%s',SYS.path,SYS.procar)
[kpnt_pos, eval, bandchar] = vasp__readPROCAR(file_procar, SYS.lorbit);

% get Fermi energy from DOSCAR
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar);
efermi = vasp__getEFermi(file_doscar);
%efermi = 0

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

% define which atoms and orbitals to plot
satoms = [5 6];
orbitals = [1:4];
magdir = 1;

% determine size of "problem"
nsatoms = length(satoms);
norbitals = length(orbitals);

% dimension of the (pdiml x pdimc) subplot matrix
% this is usually a square matrix, add one plot for the legend
pdimc = nsatoms+1; %ceil(sqrt(nsatoms))+1
pdiml = 1; %pdimc-1

for i = 1:nsatoms    
    % select the atom
    atom = satoms(i);         

    % we plot the total character of a particular atom
    bands.bchar = 0;
    for orb = 1:norbitals
        % select orbital
        orbital = orbitals(orb);
        
        % modulus of the band phases and band characters lead to
        % absolutely the same band projections; but the values of the band
        % characters are slightly bigger
        %bands.bchar = conj(bandphase(:,:,atom,orbital)).*bandphase(:,:,atom,orbital);
        bands.bchar = bandchar(:,:,atom,orbital,magdir); % wird magdir weggelassen -> magdir=1 (abwaertskompatibel !)
        

        % compile the relevant data for plotting
        bands.emax = 6; %SYS.emax;                % plotting range
        bands.emin = -14; %SYS.emin;
        bands.klabels = SYS.klabels;          % special point labels
        bands.eval = eval-efermi;         % eigen energies with Fermi level at E=0
        bands.plotbands = 1;
        bands.linespec = '-k';             % 'linespec' of the lines of the bands
        bands.kpnt_pos = kpnt_pos;        % positions of the kpoints  
        bands.charscal = .1;               % fatness of bands  = charscal * band character
        bands.charcol = color(orbital,:);

        % create the subplot
        subplot(pdiml,pdimc,i); % different subplot for each atom
        %subplot(pdiml,pdimc,orb);  % different subplot for each orbital

        % plot the bands
        vasp__plotBands(bands);
        
        % draw title
        title_str = sprintf('%s  - atom %d',SYS.ID, atom);
        title(title_str); %,'Interpreter','none');

    end
end

% add home made legend
subplot(pdiml,pdimc,nsatoms+1);        % different subplot for each atom
%subplot(pdiml,pdimc,norbitals+1);       % different subplot for each orbital
hold on
for orbital=orbitals
    ypos = 9-orbital;
    plot([0 1]',[ypos ypos]', 'Color', color(orbital,:), 'LineWidth', 3);
    text(1.15, ypos, names(orbital))
end
axis off
xlim([0 3])
ylim([-.5 8.5])

SetFontsInFigure(20)

