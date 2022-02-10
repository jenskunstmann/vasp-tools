function vasp_PlotProjBands_magdir()
% plot the projected bandstructure with respect to the magnetization
% direction

global SYS

% read PROCAR file
file_procar = sprintf('%s/%s',SYS.path,SYS.procar)
[kpnt_pos, eval, bandchar] = vasp__readPROCAR(file_procar, SYS.lorbit);

% get structure from CONTCAR
file_contcar = sprintf('%s/%s', SYS.path, SYS.contcar)
structure = cry__readCONTCAR(file_contcar);
rlv  = cry__GetRLV(structure)*2*pi;  % convert to [1/Ang]

% get Fermi energy from DOSCAR
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar);
efermi = vasp__getEFermi(file_doscar);
%efermi = 0;

% extract data
natoms = size(bandchar,3);  % number of atoms
npnts = size(kpnt_pos,1);   % number of kpoints
nbands = size(eval,2);      % number of bands
nspin = size(bandchar, 5);  % either 1 (no spin), 2(collinear), 4(non collinear)

% definition of the orbital characters:
% syntax: bandchar(kpnt,band,atom,orbital,spin)  
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
     
color_mag = [.4 .2  0; 0 0 1; 1 0 0; 0 1 0];
names_mag = {'total', '+m_x', '+m_y', '+m_z'};  

% compile the relevant data for plotting
bands.emax = SYS.emax;                % plotting range
bands.emin = SYS.emin;
bands.klabels = SYS.klabels;          % special point labels
bands.eval = eval-efermi;         % eigen energies with Fermi level at E=0
bands.plotbands = 1;
bands.linespec = '-k';             % 'linespec' of the lines of the bands
bands.kpnt_pos = kpnt_pos*rlv;        % positions of the kpoints  
bands.charscal = 0.05;               % fatness of bands  = charscal * band character

% sum up atoms and orbitals
% bandchar(kpnt,band,atom,orbital,spin)
% bchar(kpnt,band,spin)
bchar = reshape(sum(sum(bandchar,4),3), npnts, nbands, nspin);

% normalize the spin projection such that m_x^2 + m_y^2 + m_z^2 = 1
% the normalization constant is different for every state (k,n)
norm = sqrt(sum(bchar(:,:,2:4).^2, 3));
bchar(:,:,2) = bchar(:,:,2) ./ norm;
bchar(:,:,3) = bchar(:,:,3) ./ norm;
bchar(:,:,4) = bchar(:,:,4) ./ norm;


% pick numbers from the arrays
kpoint = 1
band = 17
mag_vector = reshape(bchar(kpoint, band, 2:4), 1, 3)
norm = dot(mag_vector, mag_vector)  % check norm

% uncomment to suppress the plot
% return

% select the magnetization components to p1ot
smagdir = [2 3 4];

% dimension of the (pdiml x pdimc) subplot matrix
pdiml = length(smagdir)+1; 
pdimc = 1; 

bands.bchar = 0;
for i = 1:length(smagdir)    
    
    % create the subplot
    subplot(pdiml,pdimc,i);

    % plot POSITIVE VALUES in regular color 
    bchar1 = bchar(:,:, smagdir(i));
    bchar1(bchar1 < 0) = 0;   % all negative number are set to zero
    bands.bchar = bchar1;
    bands.charcol = color_mag(smagdir(i),:);                
    vasp__plotBands(bands);

    % plot NEGATIVE VALUES in complementary color 
    bchar1 = bchar(:,:, smagdir(i));
    bchar1(bchar1 >= 0) = 0;   % all positive number are set to zero
    bands.bchar = bchar1;
    bands.charcol = [1 1 1] - color_mag(smagdir(i),:);        
    x = vasp__plotBands(bands);

    % draw title
    title_str = sprintf('%s - %s',SYS.ID, names_mag{smagdir(i)});
    title(title_str); %,'Interpreter','none');

end

% plot band numbers
%subplot(pdiml,pdimc,1);
%vasp__plotBandNumbers(bands, x, [1 0 0], [0.0 0.0])  

% add home made legend
subplot(pdiml,pdimc,length(smagdir)+1);        
hold on
for magdir = smagdir
    ypos = 4-magdir;
    plot([0 1]',[ypos ypos]', 'Color', color_mag(magdir,:), 'LineWidth', 3);
    text(1.15, ypos, names_mag(magdir))
end
axis off
xlim([0 3])
ylim([-.5 8.5])
SetFontsInFigure(12)



