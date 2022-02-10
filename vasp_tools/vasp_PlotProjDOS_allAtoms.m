function vasp_PlotProjDOS_allAtoms()
% plot projected DOS from VASP
% not working for spin-polarized calculations, so far
% 
% What this program does:
% 1. it reads in the projected dos from VASP in the following arrays
% e(pnts)                   : energy axis
% totaldos(pnts)            : total DOS
% projdos(pnts,atom,orbital): atom and orbital resolved projected DOS
%
% 2. from this 'raw' data the curves of interest are picked or constructed
% and stored in the arrays:
% pdos(pnts,curve)  : holds the data of the dos curves to be plottet
% doscolor(curve,:) : color of each curve
% doslegend{curve}  : legend of each curve
% doslstyle{curve}  : line style of each curve
% sort(curve)       : defines which curves are plottet in what order; 
%                     the order is important for the order in the legend box
%
% 3. set up the 'dos' struct and plot everythng using vasp__plotDOS()
% here plotting ranges can be adjusted


global SYS

%%% 1. reading in the projected dos from VASP 
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar);
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
atomselect = 1:natoms;

% below reshape() is repeatedly used to transform the 3-dimensional matix
% 'projdos' into the 2-dimensional one 'pdos', this is necessary

%%% first add all atoms; normalize to DOS per atom; add total DOS as 10th
%%% orbital; pdos_allat(pnts,orbital)
pdos_allat = 0;
for at=atomselect
    pdos_allat = pdos_allat + reshape(projdos(:,at,:), npnts, 9);
end
pdos_allat = pdos_allat/natoms;
pdos_allat(:,10) = totaldos;        

%%% OPTION: plot all oribtal characters 
sort = [1:9];  % the orbitals to plot  
doscolor = orbitalcolors;
doslegend = orbitalnames;
pdos = pdos_allat;

%%% OPTION: plot s + 3 p_i projected DOS
% sort = [1 4 2 3];  % the orbitals to plot  
% doscolor = orbitalcolors;
% doslegend = orbitalnames;
% pdos = pdos_allat;

%%% OPTION: SIGMA and PI states
% sort = [1 2 3];  
% doslegend = {'\sigma','\pi', 'total'};
% doscolor = [orbitalcolors(2,:);  orbitalcolors(3,:); [0 0 0]];
% %make DOS of sigma states (curve=1) 
% pdos = 0; 
% for orbital=[1 2 4]
%    pdos = pdos + pdos_allat(:,orbital);
% end
% %make DOS of pi states (curve=2)
% pdos(:,2) = pdos_allat(:,3);
% %make DOS of pi states (curve=3)
% pdos(:,3) = pdos_allat(:,10);

%%%% END: USER INPUT %%%%%
%%% now plot everything that was defined above %%%


%%% 3. set up the 'dos' struct and plot everythng using vasp__plotDOS()
dos.e      = e;
dos.ddos   = pdos;
dos.lstyle = doslstyle;     
dos.legend = doslegend;     
dos.sort   = {sort};  
dos.color  = doscolor;

dos.plotinset = 0;        % true:1: plot inset, false:0: plot only DOS

% plotting ranges
dos.emin = -18;% min(e);  %SYS.emin; %
dos.emax = 15; %max(e);  %SYS.emax; %
dos.dosmax = 0.7; %max(dos.ddos(:,1));  % take the maximum of the first DOS curve
dos.dosmin = 0;

% plotting ranges of the INSET
dos.inset_emin = -3;
dos.inset_emax = 3; 
dos.inset_dosmin = 0;
dos.inset_dosmax = dos.dosmax;

% shift inset and main plot in y direction by insetshift*dosmax
dos.insetshift = 1.5;
%%% END: user settings %%%

[mainax, insetax] = vasp__plotDOS(dos);
ylabel(mainax,'DOS (states / eV atom)');
xlabel(mainax,'E-E_F (eV)'); 

% draw title
title_str = sprintf('%s',SYS.ID);
title(mainax, title_str); %,'Interpreter','none');
SetFontsInFigure(20);

