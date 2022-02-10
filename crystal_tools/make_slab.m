function make_slab()
% atomic structure defined by
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.typeatomnum(atomID)       = atomic number of the atom 

%%% how to make a surface slab with cif2cell?
% 1. generate primitive cell
% 2. get the primitive lattice vectors a1,a2,a3
% 3. express the surface lattice vectors A1,A2,A3 in terms of a1,a2,a3
% 4. the cif2cell parameter is then --supercell=[[A1],[A2],[A3]]

% generate primitive cell
% > cif2cell Si.cif   -p vasp --vasp-format=5 
% > mv POSCAR POSCAR-primitive
% 
% generate conventional cell
% > cif2cell Si.cif   -p vasp --vasp-format=5 --no-reduce
% > mv POSCAR POSCAR-conventional
% 
% generate 100 1x1 surface slab
% > cif2cell Si.cif -p vasp --vasp-format=5 
%            --supercell=[[0,0,1],[1,-1,0],[1,1,-1]] 
%            --supercell-vacuum=[0,0,2]  
%            --supercell-postvacuum-translation=[0,0,0.5]


file_prim = '/home/jk/res/silicon/surface/initistruc/POSCAR-primitive'
file_conv = '/home/jk/res/silicon/surface/initistruc/POSCAR-conventional'
prim = cry__readCONTCAR(file_prim);
conv = cry__readCONTCAR(file_conv);

%cry__PlotAtomicStructure(crystal)
A = 5.430700;
B_prim = prim.latt/A
B_conv = conv.latt/A

% trafo
% M * B_old = B_new
% M = B_new * (B_old)^-1


% transformation primitive -> cubic conventioanl 
%B_new = [1 0 0; 0 1 0; 0 0 1]
%M = B_new * inv(B_prim)

% apply transformation as in the manual
% M = [1 1 -1; 1 -1 1; -1 1 1]
% B_new = M* B_prim

% 100 1x1 primitive surface cell
M = [0 0 1; 1 -1 0; 1 1 -1]
B_new = M* B_prim