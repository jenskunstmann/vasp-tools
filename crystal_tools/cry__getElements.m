function [atomicnum, numatoms] = cry__getElements(crystal)
% extracts data from the 'crystal' 
% ELEMENTS ARE SORTED IN ORDER OF ASCENDING ATOMIC NUMBER
%
% atomicnum(1:nelements) = atomic number of the elements in the system
% numatoms(1:nelements) = number of atoms per element

atomicnummax = 120; % there are actually only 118 elements

% sort the atomic number array in ascending order
[atomnum_sort, indexmap] = sort(crystal.atomnum, 'ascend');
% elementcount(atomicnumber) = number of atoms with that atomic number
elementcount = histc(atomnum_sort,1:atomicnummax);
% index of non-zero elements = array with the atomic numbers
atomicnum = find(elementcount); 
nelements = length(atomicnum);
%elementstr = '';
for elem =  1:nelements
    numatoms(elem) = elementcount(atomicnum(elem));    
    %elementstr = [elementstr ' ' cry__getAtomicSymbol(atomicnum(elem))];
end

