function [atomicnum, elementstr, numatoms, atompos_sort, atomnum_sort, SDmatrix_sort] = cry__GetPOSCARInfos(crystal, SDmatrix)
% extracts data from the 'crystal' structure that is required to write a
% POSCAR file
% ELEMENTS ARE SORTED IN ORDER OF ASCENDING ATOMIC NUMBER
%
% atomicnum(1:nelements) = atomic number of the elements in the system
% elementstr = string with the atomic symbols of the 'nelements' elements
% numatoms(1:nelements) = number of atoms per element
% atompos_sort(atomID, Cartesian component), sorted in ascending atomic number 
% atomnum_sort(1:natoms) = sorted atomnum array (not required for POSCAR)
% SDmatrix_sort(atomID, Cartesian component), resorted according to
%                                             the new order

atomicnummax = 120; % there are actually only 118 elements

% sort the atomic number array in ascending order
[atomnum_sort, indexmap] = sort(crystal.atomnum, 'ascend');
% elementcount(atomicnumber) = number of atoms with that atomic number
elementcount = histc(atomnum_sort,1:atomicnummax);
% index of non-zero elements = array with the atomic numbers
atomicnum = find(elementcount); 
nelements = length(atomicnum);
elementstr = '';
for elem =  1:nelements
    numatoms(elem) = elementcount(atomicnum(elem));    
    elementstr = [elementstr ' ' cry__getAtomicSymbol(atomicnum(elem))];
end

% put the atomic position in the same order as the atomic numbers
%natoms = length(crystal.atomnum);
SDmatrix_sort = [];
%for at = 1:natoms
atompos_sort = crystal.atompos(indexmap,:);    
if ~isempty(SDmatrix)
    SDmatrix_sort = SDmatrix(indexmap,:);
end
%end