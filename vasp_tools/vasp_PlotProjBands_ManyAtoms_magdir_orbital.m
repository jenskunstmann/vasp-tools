function vasp_PlotProjBands_ManyAtoms_magdir_orbital()
% plot the projected bandstructure in a matrix of subplots for many
% different atoms and show their magnetization direction

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
% syntax: bandchar(kpnt,band,atom,orbital,magdir)  
% orbital = 1 (s)   -  grey
%           2 (py)  -  red
%           3 (pz)  -  green
%           4 (px)  -  blue
%           5 (dxy) -  yellow
%           6 (dyz) -  magenta
%           7 (dz2) -  cyan
%           8 (dxz) -  orange
%           9 (dx2) -  brown
%
% magdir = 1 (total) - brown
%          2 (mx)    - blue
%          3 (my)    - red
%          4 (mz)    - green
color = [.5 .5 .5; 1 0 0; 0 1 0; 0 0 1; 1 1 0; ...
         1 0 1; 0 1 1; 1 0.5 0; .4 .2  0];
names = {'s', 'p_y', 'p_z', 'p_x', 'd_{xy}', 'd_{yz}', ...
         'd_{z^2}', 'd_{xz}', 'd_{x^2-y^2}'};  
     
color_mag = [.4 .2  0; 0 0 1; 1 0 0; 0 1 0];
names_mag = {'total', '+m_x', '+m_y', '+m_z'};  
     

% define which atoms and orbitals to plot
atom = 3;
orbitals = [1:3];
smagdir = [2 3 4];

% determine size of "problem"
%nsatoms = length(satoms);
norbitals = length(orbitals);

% dimension of the (pdiml x pdimc) subplot matrix
% this is usually a square matrix, add one plot for the legend
pdimc = 6; %ceil(sqrt(nsatoms))+1
pdiml = 1 %pdimc-1

for i = 1:norbitals    
    % select the atom
    %atom = satoms(i);         
    orbital = orbitals(i);

    % we plot the total character of a particular atom
    bands.bchar = 0;
    for magdir = smagdir
        
        % sum up all or selected orbital characters
        bchar = 0;
        %for orbital=orbitals
            bchar = bchar + bandchar(:,:,atom,orbital,magdir);
        %end                  

        % compile the relevant data for plotting
        bands.emax = 3; %SYS.emax;                % plotting range
        bands.emin = -1.5; %SYS.emin;
        bands.klabels = SYS.klabels;          % special point labels
        bands.eval = eval-efermi;         % eigen energies with Fermi level at E=0
        bands.plotbands = 1;
        bands.linespec = '-k';             % 'linespec' of the lines of the bands
        bands.kpnt_pos = kpnt_pos;        % positions of the kpoints  
        bands.charscal = 0.2;               % fatness of bands  = charscal * band character
        
        % create the subplot
        subplot(pdiml,pdimc,i); % different subplot for each atom

        % plot POSITIVE VALUES in regular color 
        bchar1 = bchar;
        bchar1(bchar1 < 0) = 0;   % all negative number are set to zero
        bands.bchar = bchar1;
        bands.charcol = color_mag(magdir,:);                
        vasp__plotBands(bands);
        
        % plot NEGATIVE VALUES in complementary color 
        bchar1 = bchar;
        bchar1(bchar1 >= 0) = 0;   % all positive number are set to zero
        bands.bchar = bchar1;
        bands.charcol = [1 1 1] - color_mag(magdir,:);        
        x = vasp__plotBands(bands);
        
        % draw title
        title_str = sprintf('%s', names{orbital});
        title(title_str); %,'Interpreter','none');

    end
end

%vasp__plotBandNumbers(bands, x, [1 0 0], [0.0 0.0])  

% add home made legend
subplot(pdiml,pdimc,norbitals+1);        % different subplot for each atom
%subplot(pdiml,pdimc,norbitals+1);       % different subplot for each orbital
hold on
for magdir = smagdir
    ypos = 4-magdir;
    plot([0 1]',[ypos ypos]', 'Color', color_mag(magdir,:), 'LineWidth', 3);
    text(1.15, ypos, names_mag(magdir))
end
axis off
xlim([0 3])
ylim([-.5 8.5])

SetFontsInFigure(20)

