function vasp_PlotProjDOS_ManyAtoms()
% plot projected DOS from VASP
% not working for spin-polarized calculations, so far
% 
% What this program does:
% see vasp_PlotProjDOS()


global SYS

%%% 1. reading in the projected dos from VASP 
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar)
[e totaldos projdos] = vasp__readprojDOSCAR(file_doscar);

% figure out the dimensions of the arrays
natoms = size(projdos,2);
npnts = size(projdos,1);


%%% 2. from the 'raw' data the curves of interest are picked or constructed 

% definition of the orbital characters:
% orbital = 1 (s)   -  grey
%           2 (py)  -  red
%           3 (pz)  -  green
%           4 (px)  -  blue
%           5 (dxy) -  yellow
%           6 (dyz) -  magenta
%           7 (dz2) -  cyan
%           8 (dxz) -  orange
%           9 (dx2) -  brown
orbitalcolors = [.5 .5 .5; 1 0 0; 0 1 0; 0 0 1; 1 1 0; ...
         1 0 1; 0 1 1; 1 0.5 0; .4 .2  0];
orbitalnames = {'s', 'p_y', 'p_z', 'p_x', 'd_{xy}', 'd_{yz}', ...
         'd_{z^2}', 'd_{xz}', 'd_{x^2-y^2}'};   
orbitallstyles = {'-','-','-','-','-','-','-','-','-'};     
doslstyle = orbitallstyles;


% define which atoms to plot
%satoms = [7,46,13,52,19,58,25,64,31,69,35,72];   % 'pureZGNR'
%satoms = [63,5,25,8,27,10,31,14,34,17,37,66];   % 'ZGNR+HT'
%satoms = [35,33,29,25,21,17,13,9,5];   % 'pureAGNR'
%satoms = [62,58,54,72,73,71,45,41,37,35];   % 'AGNR+BC' x-richtung
%satoms = [35,4,40,10,46,67,73,72,19,55,25,61];% 'ZGNR+BC' y-richtung
%satoms = [63,3,23,55,61,58,39,29,49,52,15,36];% 'ZGNR+HT+BC' y-richtung
%satoms = [48,46,59,55,57,8];% 'ZGNR+HT+BC'
%satoms = [26,3,4,27,5,30,8,31,9,34,12,35];% 'ZGNR.57' y-Richtung
%satoms = [1,26,3,25];% 'ZGNR.57' x-Richtung
%satoms = [1,26,3,25];% 'ZGNR.57+H' x-Richtung
%satoms = [1,2,3,4,5,6,7,8];% 'alpha1/res1'
%satoms = [2 5 15 24 25 26   32 33  34 35  38 83  84 85  57 58];% 'MoS2-844-Tim'
satoms = [1 5];
nsatoms = length(satoms);

% dimension of the (pdiml x pdimc) subplot matrix
% this is usually a square matrix
pdimc = 1;%ceil(sqrt(nsatoms));
pdiml = 2; %pdimc-1;

for i = 1:nsatoms

% %%% OPTION: SIGMA and PI states
% % select the atom
% atom = satoms(i); 
% sort = [1 2];  
% doslegend = {'\sigma','\pi'};
% doscolor = [orbitalcolors(2,:);  orbitalcolors(3,:) ];
% pdos_tmp = reshape(projdos(:,atom,:), npnts, 9);
% pdos = 0; 
% %make DOS of sigma states (curve=1) 
% for orbital=[1:9]
%    pdos = pdos + pdos_tmp(:,orbital);
% end
% %make DOS of pi states (curve=2) 
% pdos(:,2) = pdos_tmp(:,3);

%%% OPTION: all states
% select the atom
atom = satoms(i); 
sort = [1 2 3 ];  
nsort = 3; %length(sort);
pdos_tmp = reshape(projdos(:,atom,:), npnts, 9);
pdos = zeros(npnts, nsort);
pdos(:,1) = totaldos;
pdos(:,2) = pdos_tmp(:,1);
% sum up all p states
 for orbital=[2:4]
    pdos(:,3) = pdos(:,2) + pdos_tmp(:,orbital);
 end
% sum up all d states
 for orbital=[5:9]
    pdos(:,4) = pdos(:,3) + pdos_tmp(:,orbital);
 end
doslegend = {'total','s', 'p', 'd'}; 
doscolor = [0 0 0; 1 0 0; 0 0 1; 0 1 0; ...
         1 0 1; 0 1 1; 1 0.5 0; .4 .2  0];
% %make DOS of pi states (curve=2) 
% pdos(:,2) = pdos_tmp(:,3);


%%% 3. set up the 'dos' struct and plot everythng using vasp__plotDOS()
dos.e      = e;
dos.ddos   = pdos;
dos.lstyle = doslstyle;     
dos.legend = doslegend;     
dos.sort   = {sort};  
dos.color  = doscolor;

dos.plotinset = 0;        % never plot inset, does not work in a series of subplots

% plotting ranges
dos.emin =  -14; %min(e); % SYS.emin;
dos.emax =  10;%max(e); % SYS.emax;
dos.dosmax = 1;%max(dos.ddos(:,1));  % take the maximum of the first DOS curve
dos.dosmin = 0;

% plotting ranges of the INSET
dos.plotinset = 0;
dos.inset_emin = -2;
dos.inset_emax = 3; 
dos.inset_dosmin = 0;
dos.inset_dosmax = dos.dosmax;

% shift inset and main plot in y direction by insetshift*dosmax
dos.insetshift = 1.0;
%%% END: user settings %%%

% create the subplot
subplot(pdiml,pdimc,i);
[mainax, insetax] = vasp__plotDOS(dos);

% draw title
title_str = sprintf('%s  - atom %d',SYS.ID, atom);
title(mainax, title_str); %,'Interpreter','none');
ylabel(mainax,'DOS (states / eV atom)');
xlabel(mainax,'E-E_F (eV)');

%SetFontsInFigure(20);

end
