function alltogether = cry__MergeAtomicStructure(varargin)
% merge atomic structures into one
% the lattice vectors of the final structure is the one of the first
% structure
%
% USAGE: alltogether = MergeAtomicStructure(Patch_1, Patch_2, ...)
%
% atomic structure defined by
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.atomnum(atomID)       = atomic number of the atom 

% parse input arguments, nargin includes the explicitly specified arguments
if (nargin < 2)
    error('Too few input arguments');
end
npatches = nargin;
for ipat = 1:npatches    
    patch(ipat) = varargin{ipat};
    numat(ipat) = size(patch(ipat).atompos,1);
end

% merge the atomic structures
at = 0;
for ipat = 1:npatches    
    % atomic positions
    at_start = at + 1;
    at_end   = at + numat(ipat);
    alltogether.atompos(at_start:at_end,:) = patch(ipat).atompos;
    alltogether.atomnum(at_start:at_end)   = patch(ipat).atomnum;
    at = at_end;     
end

% lattice vectors = the one of the first structure
alltogether.origin = patch(1).origin;
alltogether.latt = patch(1).latt;

end

