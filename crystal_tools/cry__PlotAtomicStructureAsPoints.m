function cry__PlotAtomicStructureAsPoints(crystal, varargin)
% plots the atomic structure and the unit cell from a crystal data
% structure,  not as spheres but as symbols
% their atomic numbers; cell repetitions are possible, different species
% have different colors and radii, atomic position are displayed in Angstrom
% the supercell is defined by the optional arguments cmin, cmax
% color and radius can also be defined differently for each atom ID
%
% usage: PlotAtomicStructure(crystal, [cmin, cmax, [atomicradius, atomiccolor], ['AtomIDBalls']])
% if the optional arguments are empty, i.e. [], the default values are used      
%
% atomic structure defined by
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.atomnum(atomID)       = atomic number of the atom 
%
% the atomic numbers can also be determined in xcrysden, 
% but ONLY IF one unit cell is drawn
% USES: the function arrow3()

%%% default plotting options
% as these options are only relevant for this routine it is not necessary
% to create a global atomic properties structure
atomicnummax = 120; % there are actually only 118 elements
% some generic options
atomicradius = linspace(0.1, 0.5, atomicnummax); % radius for each atom type
atomiccolor  = [linspace(.5, 0, atomicnummax)' ...
                linspace(.5, 0, atomicnummax)' ...
                linspace(.5, 0, atomicnummax)'];   % color string for each atom type
% override the default parameters for special cases; to be extended if
% needed
atomicradius(1)  = .3; atomiccolor(1,:)  = [0 0 0]; % H
atomicradius(5)  = .3; atomiccolor(5,:)  = [0.3 0.3 1]; % B
atomicradius(6)  = .5; atomiccolor(6,:)  = [0 0 0]; % C
atomicradius(16) = .4; atomiccolor(16,:) = [1 1 0]; % S
atomicradius(42) = .6; atomiccolor(42,:) = [0 0 1]; % Mo         

% plot atom numbers or not
plot_atomID = false;

            
% position of the atomic number text relative to the atomic position and relative to the atomic radius away from the atomic sphere  
shift = -1.05*[1 1 1];               
% supercells, range of plotting
cmin = [0 0 0];     
cmax = [0 0 0];

% parse input arguments to modify plotting options
% only overwrite if non-empty value is given
AtomIDBalls = false;
switch(nargin)
    case(1)

    case(3)
        if ~isempty(varargin{1})
            cmin = varargin{1};
        end        
        if ~isempty(varargin{2})
            cmax = varargin{2};
        end          
    case(5)
        if ~isempty(varargin{1})
            cmin = varargin{1};
        end        
        if ~isempty(varargin{2})
            cmax = varargin{2};
        end   
        if ~isempty(varargin{3})
            atomicradius = varargin{3};
        end
        if ~isempty(varargin{4})
            atomiccolor = varargin{4};
        end        
    case(6)
        % atomicradius and atomiccolor are not given for every atom TYPE but
        % for every atomID
        AtomIDBalls = true;
        if ~isempty(varargin{1})
            cmin = varargin{1};
        end        
        if ~isempty(varargin{2})
            cmax = varargin{2};
        end   
        if ~isempty(varargin{3})
            atomicradius = varargin{3};
        end
        if ~isempty(varargin{4})
            atomiccolor = varargin{4};
        end        
    otherwise
        error('wrong number of input arguments')
end


% find out number of atoms
natoms = size(crystal.atompos,1);

% make a sphere surface
[X0,Y0,Z0] = sphere(6); 

% add the spheres to the plot
hold on
for i1 = cmin(1):cmax(1)     % loop over the supercell lattice vectors
    for i2 = cmin(2):cmax(2)
        for i3 = cmin(3):cmax(3)            
            for at = 1:natoms       % loop over the atoms
                    
                % for readability of the code we define a few variables              
                pos = [i1 i2 i3]*crystal.latt + crystal.atompos(at,:);
                atomnum = crystal.atomnum(at);
                 
                % ser color and radius depending on atom type OR atom ID
                if AtomIDBalls
                    color = atomiccolor(at,:);
                    rad   = atomicradius(at);  
                else                    
                    color = atomiccolor(atomnum,:);
                    rad   = atomicradius(atomnum);  
                end
                shi = shift*rad; 
                    
                % put the atomic spheres and atomic numbers at their place                
%                 surf(rad*X0 + pos(1), rad*Y0 + pos(2), rad*Z0 + pos(3), ...
%                     'FaceColor',color, 'EdgeColor','none','AmbientStrength',.5);

                plot3(pos(1), pos(2), pos(3), 'ok', 'Color', color, 'MarkerSize', rad)
                
                if plot_atomID
                    text(pos(1)-shi(1), pos(2)+shi(2), pos(3)+shi(3), mat2str(at))%, ...
                         %'FontSize',12); %, 'FontName','Arial');
                end
                    
            end;
        end;
    end;
end;


plotcells(crystal.latt, crystal.origin, cmin, cmax, 'black')
% 
% % display options
daspect([1 1 1])    % aspect ratio
%axes3;          % arrows for x,y,z axes
% 
% DEFINE VIEW
%view(-16,23);    % viewpoint: azimuth = 30 and elevation = 40
%camzoom(1.3);    % zoom in
%camproj perspective;    % perspective projection

% DEFINE LIGHTING OPTIONS 
% camlight(az,el): az=phi:-180..180  el=theta:-90(bottom)..90(top)
% camlight(-130,-60);    % creates two lights; one from above
% camlight(50,60)   % the other from below
% lighting phong;  % lighting mode: phong, none, flat, gouraud
% 
% shg;            % bring to front
% rotate3d;       % turn on rotation mode


