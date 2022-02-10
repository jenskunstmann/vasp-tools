function vasp_PlotProjDOS()
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

e_Fermi = 0.5; 

% figure out the dimensions of the arrays
natoms = size(projdos,2);
npnts = size(projdos,1);


%%% 2. from the 'raw' data the curves of interest are picked or constructed 

% definition of the orbital characters:
% orbital = 1 (s)   -  grey
%           2 (py)  -  red
%           3 (pz)  -  blue
%           4 (px)  -  green
%           5 (dxy) -  yellow
%           6 (dyz) -  magenta
%           7 (dz2) -  cyan
%           8 (dxz) -  orange
%           9 (dx2) -  brown
orbitalcolors = [.5 .5 .5; 1 0 0; 0 0 1; 0 1 0; 1 1 0; ...
         1 0 1; 0 1 1; 1 0.5 0; .4 .2  0];
orbitalnames = {'s', 'p_y', 'p_z', 'p_x', 'd_{xy}', 'd_{yz}', ...
         'd_{z^2}', 'd_{xz}', 'd_{x^2-y^2}'};   
orbitallstyles = {'-','-','-','-','-','-','-','-','-'};     
doslstyle = orbitallstyles;


% below reshape() is repeatedly used to transform the 3-dimensional matix
% 'projdos' into the 2-dimensional one 'pdos', this is necessary

%%% OPTION: plot s + 3 p_i projected DOS of a specific atom
% atom = 1;          % select the atom to plot:
% sort = [1 4 2 3];  % the orbitals to plot  
% doscolor = orbitalcolors;
% doslegend = orbitalnames;
% pdos = reshape(projdos(:,atom,:), npnts, 9);

% integrated characters:

%%% OPTION: orbital projected dos, integrated over ALL atoms
% = sum up all atomic contributions   
% sort = [1:4];    % the orbitals to plot  
% doscolor = orbitalcolors;
% doslegend = orbitalnames;
% pdos = 0; 
% for atom = 1:natoms    
%     pdos = pdos + reshape(projdos(:,atom,:), npnts, 9)/natoms; 
% end

%%% OPTION: SIGMA and PI states
% atom = 2;
% sort = [1 2];  
% doslegend = {'\sigma','\pi'};
% doscolor = [orbitalcolors(2,:);  orbitalcolors(3,:) ];
% %pdos_tmp = reshape(projdos(:,atom,:), npnts, 9);
% pdos_tmp = pdos; % copy integrated pdos from above
% pdos = 0; 
% %make DOS of sigma states (curve=1) 
% for orbital=[1 2 4]
%    pdos = pdos + pdos_tmp(:,orbital);
% end
% %make DOS of pi states (curve=2) 
% pdos(:,2) = pdos_tmp(:,3);

%%% OPTION: self-made state combinations
%atom = 8;
%sort = [1 2 3];  
%doslegend = {'s + p_z', 'p_x', 'p_y'};
%doscolor = [orbitalcolors(3,:);  orbitalcolors(4,:); orbitalcolors(2,:) ];
%pdos_tmp = reshape(projdos(:,atom,:), npnts, 9);
%pdos = 0; 
% make DOS of sigma states (curve=1) 
%for orbital=[1 3]
%    pdos = pdos + pdos_tmp(:,orbital);
%end
% make DOS of pi states (curve=2) 
%pdos(:,2) = pdos_tmp(:,4);
%pdos(:,3) = pdos_tmp(:,2);


%%% OPTION: total dos of different atoms
%= sum up all orbital characters
%which atoms to plot
pdos = 0;
for orbital=1:4
   pdos = pdos + reshape(projdos(:,:,orbital), npnts, natoms); 
end

pdos_select = zeros(npnts, 3);
% sum up S
for atom=[1:20]
   pdos_select(:,1) = pdos_select(:,1) + pdos(:,atom); %/length(atom);
end

% sum up Sn
for atom=[21:36]
   pdos_select(:,2) = pdos_select(:,2) + pdos(:,atom); %/length(atom);
end

% sum up Pb
for atom=[37:40]
   pdos_select(:,3) = pdos_select(:,3) + pdos(:,atom); %/length(atom);
end

% projdos(pnt,atom,orbital)
%pdos_select(:,1) = projdos(:,1,1); 
%pdos_select(:,2) = projdos(:,3,1); 


sort = [1 2 3]; 
% here 'pdos' has a large dimension: pdos(1..npnts,1..natoms)
% therefore we have to provide doscolor, doslstyle and doslegend with the
% same large dimension, in this loop we just define the elements that are
% needed
for i = sort
    %doscolor(i,:) = [i/natoms 1-i/natoms rand];
    doslstyle{i} = '-';
    %doslegend{i} = sprintf('total of atom %d', i);
end
doscolor(sort,:) = orbitalcolors([2 3 9],:);
doslegend = {'S', 'Sn', 'Pb'};


%%% OPTION: total dos of the sum of different atoms
% = sum up all orbital characters and the specified atoms
% which atoms to plot
% sort = [1]; 
% atom = 1; % no function here, just to make the title work
% atnum = [2 5 15 24 25 26   32 33  34 35  38 83  84 85  57 58];% 'MoS2-844-Tim'
% nprojatoms = length(atnum);
% pdos_1 = 0;
% for orbital=1:9
%    pdos_1 = pdos_1 + reshape(projdos(:,:,orbital), npnts, natoms)/9; 
% end
% % here 'pdos_1' has a large dimension: pdos_1(1..npnts,1..natoms)
% 
% pdos = 0;
% for i=atnum
%     pdos = pdos + pdos_1(:,i)/nprojatoms;
% end
% doscolor = [0 0 0];
% %doslstyle = '-';
% doslegend = {'total of GB'};


%%% OPTION: specific orbital projected dos of different atoms
%orbital=3;          % select the orbital to plot
%sort = [1:natoms];  % list of atoms to plot
%doscolor = orbitalcolors;
%pdos = reshape(projdos(:,:,orbital), npnts, natoms); 
% make a legend for each element of 'sort'
%for i = 1:size(sort,2)
%    doslegend{i} = sprintf('%s of atom %d', orbitalnames{orbital}, sort(i));
%end

%%%% END: USER INPUT %%%%%
%%% now plot everything that was defined above %%%


%%% 3. set up the 'dos' struct and plot everythng using vasp__plotDOS()
dos.e      = e-e_Fermi;
dos.ddos   = pdos_select; %pdos;
dos.lstyle = doslstyle;     
dos.legend = doslegend;     
dos.sort   = {sort};  
dos.color  = doscolor;

dos.plotinset = 0;        % true:1: plot inset, false:0: plot only DOS

% plotting ranges
dos.emin =  -3; %SYS.emin;
dos.emax =  3; %SYS.emax;
dos.dosmax = 20; %max(dos.ddos(:,1));  % take the maximum of the first DOS curve
dos.dosmin = 0;

% plotting ranges of the INSET
dos.inset_emin = -3;
dos.inset_emax = 3; 
dos.inset_dosmin = 0;
dos.inset_dosmax = dos.dosmax;

% shift inset and main plot in y direction by insetshift*dosmax
dos.insetshift = 1.0;
%%% END: user settings %%%

[mainax, insetax] = vasp__plotDOS(dos);
%ylabel(mainax,'DOS (states / eV)');
ylabel(mainax,'DOS (a.u.)');
xlabel(mainax,'E-E_F (eV)'); 

% draw title
%title_str = sprintf('%s - atom %d, %s character',SYS.ID, atom, orbitalnames{orbital});
%title_str = sprintf('%s - atom %d',SYS.ID, atom);
%title(title_str); %,'Interpreter','none');

SetFontsInFigure(20);

