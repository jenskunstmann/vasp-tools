function vasp_Plot2Dbands()
% make a 2D band structure plot
% 
% KPOINTS file has to be generated first with vasp_generate_KPOINTS()
% because VASP uses shifted grids 
% VASP options: setup for DOS caclculation, ISYM = -1, ISMEAR = 0 (no
% tetahedron method, because we do not provide the tetrahedra in the
% KPOINTS)  
% default VASP setting: NKDIMD = 10000 = max. # of irred. k-points 

global SYS

% user input
bandnrs = [10];  % band numbers to plot
emin= -2;        % energy range
emax = 0;

file_eigenval = sprintf('%s/%s', SYS.path, SYS.eigenval_dos)
[kpnt_pos eval_up eval_down] = vasp__getEigenvalues(file_eigenval);

% kpnt_pos(kpoint, component): position of the k-point in reduced
%                              coordinates,  i.e. in units of the 
%                              reciprocal lattice vectors
% eval_up/down(kpoint, band) : band energies at a particular k-point 

% get Fermi energy from DOSCAR
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar);
efermi = vasp__getEFermi(file_doscar)

% get reciprocal lattice vectors (rlv) from atomic structure from CONTCAR
% [rlv] = [1/Ang]
file_contcar = sprintf('%s/%s', SYS.path, SYS.contcar)
crystal = cry__readCONTCAR(file_contcar);
rlv  = cry__GetRLV(crystal)*2*pi;  % convert to [1/Ang]

% Cartesian coordinates
kpnt_pos = kpnt_pos*rlv; % Cartesian positions of the kpoints in units of [1/Ang]

% substract fermi energy
if(isempty(eval_down))     % non-spin-polarized if 'eval_down' is empty
    %display('non-SP')
    eval = eval_up - efermi;     
else                             % spin-polarized
    % spin-up bands
    eval = eval_up - efermi;     
    % spin-down bands
    %bands.eval = eval_down - efermi;   
end;

% extract data
nkpnts = size(kpnt_pos,1)   % number of kpoints
nbands = size(eval,2)    % number of bands
kmin1 = min(kpnt_pos(:,1));
kmax1 = max(kpnt_pos(:,1));
kmin2 = min(kpnt_pos(:,2));
kmax2 = max(kpnt_pos(:,2));


% surfaces in Matlab are represented by three 2D matrices: X,Y,Z
% (i,j) -> X
% (i,j) -> Y
% (i,j) -> Z
% representing the three Cartesian coordinates of a parametric surface

npoints_int = 300; % number of points per direction for interpolation
kx = linspace(kmin1,kmax1,npoints_int);
ky = linspace(kmin2,kmax2,npoints_int);
[X,Y] = meshgrid(kx,ky);

% start plotting
hold on
%colormap summer;
colormap Jet;
%color = ['g' 'b' 'r'];


for i = bandnrs    
    % use the original data, only works in fractional coordinates
    %Z = reshape(eval(:,i), 99, 99);
    
    % interpolate the original data points
    % the interpolation also creates a uniform Cartesian mesh, i.e. the
    % new mesh is different from the old mesh
    Z = griddata(kpnt_pos(:,1), kpnt_pos(:,2), eval(:,i), X,Y,'linear');
    %size(Z)    
    
    % plot surface
    surf(X,Y,Z, ...        
         'EdgeColor','none','FaceColor','interp','BackFaceLighting','lit',...
         'AmbientStrength',0.5, 'FaceAlpha', 1);
    %colormap Gray;
    
    % plot contour
    contour3(X,Y,Z,15,'-w');
    
    % superimpose the original data points 
    %plot3(kpnt_pos(:,1), kpnt_pos(:,2), eval(:,i),'.k');        
end

% add reciprocal unit cell
newrlv = rlv;
newrlv(3,3) = abs(emax-emin);
plotcells(newrlv, [0 0 emin], [0 0 0], [0 0 0], 'black')


% end plotting
hold off


% set figure options
% Further plotting options can be set in the 
% Menue in the Figure window: File -> Export Setup
colorbar
daspect([1 1 1])
xlabel('k_x [1/Ang]');
ylabel('k_y [1/Ang]');
zlabel('E [eV]');
lighting phong;     % lighting mode: phong, none, flat, gouraud
camlight(-130,-60); % creates two lights; one from above
camlight(50,60)     % the other from below the plane
view(0, 90); %view(-35,30)        % camera position
shg;                % bring to front
box on             % axes boxed
rotate3d;           % turn on mouse rotation


% save band energies in an ASCII file, can be read in xmgrace as block data
format = 'kpnt_pos(kpointID, Cartesian component) [1/Ang], eval(kpointID, band number) [eV]'
%save('bands.mat', 'kpnt_pos', 'eval', 'format')
