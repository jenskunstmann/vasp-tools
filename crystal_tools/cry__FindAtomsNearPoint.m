function [atoms, dist_array, dist_vecs] = cry__FindAtomsNearPoint(crystal, points, distance, varargin)
% find all atoms of the crystal structure that are within 'distance'
% from the set of points, search dimension can be specified
%
% USAGE: cry__FindAtomsNearPoint(crystal, points, distance [, dim])
%
% crystal  = 'crystal' structure
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.atomnum(atomID)       = atomic number of the atom 
%
% points(numvec,comp) = Cartesian coordinates of the set of points ; 
%                        numvec = 1,..,ncorners; 
%                        comp=1,2,3, Cartesian component
% dim() = array defining the dimension(s) to consider, default = 1:3, 
%
% atoms(:)      = vector with atom IDs that are found 
% dist_array(:) = corresponding distances

% parse input
switch(nargin)
    case(3)
        dim = 1:3;      % make a 3D search
        
    case(4)
        dim = varargin{1}; % variable dimension(s)
        
    otherwise
        error('Wrong number of arguments')        
end

% N = number points
[N, three] = size(points);
natoms = length(crystal.atomnum);

if N <1
    error('Need at least one point.')
end

% % construct the edge vectors connecting the corner points
% for c = 1:(N-1)
%    edge(c,:) = points(c+1,:) - points(c,:);  
% end
% edge(N,:) = points(1,:) - points(N,:);

% find the atoms, 
atoms = 0;
foundat = 0;
for at = 1:natoms
    
    for pnt = 1:N  % loop over (corner points, edge) pairs
        % create short variables for readability
        distvec = crystal.atompos(at,dim) - points(pnt,dim);
        dist    = norm(distvec);                
        
        if dist <= distance
            foundat = foundat + 1;
            atoms(foundat) = at; 
            dist_array(foundat) = dist;
            dist_vecs(foundat,:) = distvec;
            break % edge loop and consider next atom
        end        

    end
    
end      
    
end
  
