function vasp_PlotProjBands_2layers()
% projected band structures for bilayers and heterostructures
% all atoms on the individual layers are summed up, orbital projections
% are still possible

global SYS

% read PROCAR file
file_procar = sprintf('%s/%s',SYS.path,SYS.procar)
[kpnt_pos, eval, bandchar] = vasp__readPROCAR(file_procar, SYS.lorbit);
% or read in before and get it from workspace
%global kpnt_pos eval bandchar

% get Fermi energy from DOSCAR
efermi_shift = 0.0;
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar)
efermi = vasp__getEFermi(file_doscar) + efermi_shift;

% get structure from CONTCAR
file_contcar = sprintf('%s/%s', SYS.path, SYS.contcar)
structure = cry__readCONTCAR(file_contcar);
rlv  = cry__GetRLV(structure)*2*pi;  % convert to [1/Ang]

% extract data
natoms = size(bandchar,3);  % number of atoms
npnts = size(kpnt_pos,1);   % number of kpoints
nbands = size(eval,2);    % number of bands

% determine center of mass in z-direction = point to split the layers
mean_z = mean(structure.atompos(:,3))

%%%% split structure in upper and lower layer
[~, ~, atomID_up, atomID_down] = hs__SplitBilayer(structure, mean_z);
natoms_up   = length(atomID_up)
natoms_down = length(atomID_down)

Se_atomID = cry__FindAtomsByType(structure, 34);
[UpperSeAtomIDs, LowerSeATomIDs] = FindUpperAndLower(structure, Se_atomID);



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
bands.linespec = '-k';             % 'linespec' of the lines of the bands
bands.klabels = SYS.klabels;          % special point labels
%bands.kpnt_pos = kpnt_pos;        % positions of the kpoints
bands.kpnt_pos = kpnt_pos*rlv;    % Cartesian positions of the kpoints in units of [1/Ang]
bands.eval = eval-efermi;         % eigen energies with Fermi level at E=0
bands.plotbands = 1;    
%bands.color = [1 1 1]; 
bands.charscal = 0.1;               % fatness of bands  = charscal * band character  

% list of oribitals to be summed up for the plot
sumorbitals = [1:9];

SYS.emin = -3;
SYS.emax = 3;

%%%% upper layer
subplot(1,2,1)
bands.bchar = 0;
for atom = atomID_up  % UpperSeAtomIDs
    for orbital=sumorbitals
        bands.bchar = bands.bchar + bandchar(:,:,atom,orbital);
    end    
end

% interlayer hybridization
kpoint_nr = 1
band_nr = 18
totalcharacter_up = bands.bchar(kpoint_nr,band_nr)

bands.charcol = [0 100/256 0];
bands.emin = SYS.emin;
%bands.emin = 1.2; 
bands.emax = SYS.emax;                 % plotting range
%bands.emax = 1.65;                 % plotting range
vasp__plotBands(bands);

% draw title
%title_str = sprintf('%s, outer Se atoms',SYS.ID);
title_str = sprintf('%s, WSe2',SYS.ID);
title(title_str); %,'Interpreter','none');

%%%% lower layer
subplot(1,2,2)
bands.bchar = 0;
for atom = atomID_down  % LowerSeATomIDs
    for orbital=sumorbitals
        bands.bchar = bands.bchar + bandchar(:,:,atom,orbital);
    end    
end
bands.charcol = [0 0 1]; %[1 0.6 0];
bands.emin = SYS.emin;
%bands.emin = 1.2; %-1.35; 
bands.emax = SYS.emax;                
%bands.emax = 1.65; %-0.9;                % plotting range
virtax = vasp__plotBands(bands);  % plot bands and return the virtual axis
ylabel('')

vasp__plotBandNumbers(bands, virtax, [1 0 0], [0.0 0.0])    
band_nr = 18;
% nr_of_points4fit = 5;
% m = el__getBandMass(virtax, bands, band_nr, nr_of_points4fit)

% interlayer hybridization
totalcharacter_down = bands.bchar(kpoint_nr,band_nr)

% draw title
%title_str = sprintf('%s, inner Se atoms',SYS.ID);
title_str = sprintf('%s, MoS2',SYS.ID);
title(title_str); %,'Interpreter','none');

SetFontsInFigure(20);


function [UpperAtomIDs, LowerATomIDs] = FindUpperAndLower(struc, SelectedAtomIDs)
% find the atomIDs of upper and lower atoms wrt the mean z position;
% seach only selected atoms, defined by 'SelectedAtomIDs'

mean_z = mean(struc.atompos(SelectedAtomIDs,3));


UpperAtomIDs = []; uppercnt = 0;
LowerATomIDs = []; lowercnt = 0;
for at = SelectedAtomIDs
    if(struc.atompos(at,3) > mean_z)
        uppercnt = uppercnt + 1;
        UpperAtomIDs(uppercnt) = at;
    else
        lowercnt = lowercnt + 1;
        LowerATomIDs(lowercnt) = at;        
    end    
end

