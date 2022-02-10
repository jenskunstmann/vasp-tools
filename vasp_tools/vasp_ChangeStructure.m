function vasp_ChangeStructure()
% read POSCAR/CONTCAR file modify and save again
%
% 'crystal' datastructure
% atomic structure defined by
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.atomnum(atomID)       = atomic number of the atom 
%

global SYS

%%%%% read DOSCAR file
%file_contcar = '/home/jens/prog/python/test/sheet_eps/CONTCAR-m0.02'
file_contcar = sprintf('%s/%s', SYS.path, SYS.contcar)

outfile = sprintf('%s_mod',file_contcar);
crystal = cry__readCONTCAR(file_contcar);

%%%%% print out details
natoms = length(crystal.atomnum)

%%%% general data

% dimensions
file_contcar
[A,B,C] = getLattConst(crystal.latt)
% angle definitions follow chrystallographic conventions
alpha = vecangle(crystal.latt(2,:),crystal.latt(3,:))
beta  = vecangle(crystal.latt(1,:),crystal.latt(3,:))
gamma = vecangle(crystal.latt(1,:),crystal.latt(2,:))

% % finding the primitive cell from a supercell 
% T = [1 1  0; -1 1 0; 0 0 1]         % supercell transformation matrix (as in cif2cell)
% primcell_back = crystal.latt/T      % inverted supercell trafo
% [A,B,C] = getLattConst(primcell_back)


% 3x3x1 supercell
%cry_out = cry__MakeSuperCell(crystal, [0 0 0], [2 2 0]);
%cry_out = cry__RemoveAtoms(cry_out, [1])

%crystal.atomnum


% % this is the algorithm for automatic k-point generation in VASP
% rlv = makeRLV(crystal.latt)
% rlv = cry__GetRLV(crystal)
% 
% l=24
% l*norm(rlv(1,:))
% 
% N1 = floor(max(1, l*norm(rlv(1,:)) + 0.5))
% N2 = floor(max(1, l*norm(rlv(2,:)) + 0.5))
% N3 = floor(max(1, l*norm(rlv(3,:)) + 0.5))
% 
% rb1 = [0:(N1-1)]/N1
% rb2 = [0:(N2-1)]/N2
% rb1 = [0:(N3-1)]/N3

%%%%% change things

% remove atoms
% atomIDs = [24 27 5 11 49 53 29 33 4 10 41 46];
% crystal = cry__RemoveAtoms(crystal, atomIDs);


% change hexagonal cell
% A = A - 6;
% latt2d = A*[sqrt(3)/2  -.5; sqrt(3)/2 .5];
% crystal.latt(1:2,1:2) = latt2d;
% crystal.latt

% apply in plane strain
% strain = +10 % in % 
% 
% % this code is stupid because the atoms are not scaled
% crystal.latt  % original lattice
% new_latt2d = crystal.latt(1:2,1:2)*(1+strain/100) % scale 2D vectors
% crystal.latt(1:2,1:2) = new_latt2d; % apply changes
% crystal.latt  % new lattice vectors


%crystal.latt(1:2,1:2) = latt2d;


% in relative coordinates
%crystal.atompos/crystal.latt
%crystal.atompos

%interlayer separation in 2H MoS2
%ILSep_old = crystal.atompos(6,3) - crystal.atompos(5,3)
%estimated increase of IL separation for 30 deg. twisted 2L
%DeltaSS = -0.16;

% new inter-layer separation
%ILSep_new = ILSep_old + DeltaSS

% %shift parts of the atoms that are above a certain z plane
% cuttingline_z = 6.3; % Ang
% shift = [0 0 DeltaSS];
% for at = 1:natoms
%     if crystal.atompos(at,3) > cuttingline_z
%         crystal.atompos(at,:) = crystal.atompos(at,:) + shift;
%     end
% end
% 
% % check interlayer separation
% ILSep_actual = crystal.atompos(6,3) - crystal.atompos(5,3)

% shift parts of the atoms
% for at=1:natoms
%     % shift atoms in the upper half (with respect to z) to the lower half
%     if structure.atompos(at,3) > 0.5
%         structure.atompos(at,3) = structure.atompos(at,3) - 1;
%     end
%     % shift right parts to the left
%     if structure.atompos(at,1) > 0.6
%         structure.atompos(at,1) =  structure.atompos(at,1) -1;
%     end
% end


% apply a global shift to the Cartesian coordinates
% shift = [0, 0, .45]*crystal.latt; % shift in units of the lattice vectors
% for at = 1:natoms
%     crystal.atompos(at,:) = crystal.atompos(at,:) + shift;
%     %cartpos(at,:) = cartpos(at,:) + shift;
% end;


% change lattice constants
% DeltaLatt = [0 0 0; 0 0 0; 0 0 11.5];
% crystal.latt = crystal.latt + DeltaLatt

% % set C cell dimentsion to a specific value
% Cold = crystal.latt(3,:);
% Cnew = [0 0 14];
% crystal.latt(3,:) = Cnew;
 
% % apply a shift again
% % [A,B,C] = getLattConst(crystal.latt);
% a = cartpos(25,:)
% b = crystal.latt(1,:)
% c = [0 0 7]
% d = [-1 6 0]
% shift = -(a-b)+c+d
% for at=1:natoms
%     cartpos(at,:) = cartpos(at,:) + shift;
% end


% structure.format = 'Cartesian';
% structure.atompos = cartpos;

%multiply(structure, [1 2 1])

% find angle between twisted bilayers
% % atnum_top = [16 17];
% % atnum_bottom = [25 107];
% % 
% % vec_top = crystal.atompos(atnum_top(2),:) - crystal.atompos(atnum_top(1),:)


%%%% display result
%cry__PlotAtomicStructure(crystal);

%%%% save in new file
%cry__WritePOSCAR(outfile, crystal)

end


function [A,B,C] = getLattConst(latt)
% determine lattice constants from the lattice vectors
A = norm(latt(1,:),2);
B = norm(latt(2,:),2);
C = norm(latt(3,:),2);

end


