function atoms = cry__FindAtomsInPolygon(crystal, corners)  
% find all atoms of the crystal structure that are inside the
% polygon; the corner points must be given in CONTER-CLOCKWISE order,
% closing one loop when going from the first to the last point
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

if N <3
    error('Polygon needs at least three corners.')
end

% construct edge vectors connecting the corner points, they should form a
% closed counter-clockwise loop 
for c = 1:(N-1)
   edge(c,:) = corners(c+1,:) - corners(c,:);
end
edge(N,:) = corners(1,:) - corners(N,:);

% find the atoms, only works for atom inclusion, changing the order of the
% cornerpoints to clockwise will NOT invert the search
atoms = 0;
foundat = 0;
for at = 1:natoms
    atominside = true; % reset label
    
    for ed = 1:N  % loop over (corner points, edge) pairs
        % consider atomic position relative to the corner point defining the
        % origin of the current edge vector
        pos = crystal.atompos(at,:) - corners(ed,:);
        
        % negative z component of the cross product means that the atom
        % lies RIGHT from the current edge vector edge(i,:), i.e. outside
        % of the polygon
        cp = cross(edge(ed,:),  pos);
        if cp(3) < 0
            atominside = false;
            %break % the loop over edges
        end
    end
    
    if atominside        
        foundat = foundat + 1;
        atoms(foundat) = at;
    end        
    
end

end
