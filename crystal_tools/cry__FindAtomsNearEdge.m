function atoms = cry__FindAtomsNearEdge(crystal, corners, distance)
% find all atoms of the crystal structure that are within 'distance'
% from the edges of the polygon
% 
% crystal  = 'crystal' structure
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.atomnum(atomID)       = atomic number of the atom 
%
% corners(numvec,comp) = corner points of polygon; 
%                        numvec = 1,..,ncorners; 
%                        comp=1,2,3, Cartesian component
%
% atoms() = vector with atom IDs that are found 

% N = number of corners of polygon
[N, three] = size(corners);
natoms = length(crystal.atomnum);

if N <2
    error('Need at least two corner points')
end

% construct the edge vectors connecting the corner points
for c = 1:(N-1)
   edge(c,:) = corners(c+1,:) - corners(c,:);  
end
edge(N,:) = corners(1,:) - corners(N,:);

% find the atoms, 
atoms = 0;
foundat = 0;
for at = 1:natoms
    
    for ed = 1:N  % loop over (corner points, edge) pairs
        
        % create short variables for readability
        p = crystal.atompos(at,:);        
        a = corners(ed,:);
        n = edge(ed,:)/norm(edge(ed,:));
        
        % distance of atom from edge
        dist = norm((a-p) - (dot((a-p), n) * n));
        
        % debugging output
        %[at ed dist]
        
        if dist <= distance
            foundat = foundat + 1;
            atoms(foundat) = at; 
            break % edge loop and consider next atom
        end        

    end
    
end      
    
end
  
