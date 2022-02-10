function atomicsymbol = cry__getAtomicSymbol(atomicnumber)
%  return the atomic symbol for the given atomic number

% get atomsymbols to atomnumber mapping
atom2num = cry__getAtom2Num();
atom2num_len = size(atom2num,1);

% search for the cell array index that gives us the atomic symbol 
found = false;
for index = 1:atom2num_len
    if  atomicnumber == cell2mat(atom2num(index,2))
        found = true;
        atomicsymbol = cell2mat(atom2num(index,1));
        break;
    end
end

if found == false
    error(sprintf('Atomic symbol for atomic number %d not found.\nCheck routine cry__getAtom2Num().', atomicnumber))
end