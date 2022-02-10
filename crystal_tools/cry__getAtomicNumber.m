function atomicnumber = cry__getAtomicNumber(atomicsymbol)
% return the atomic number for the given atomic symbol

% get atomsymbols to atomnumber mapping
atom2num = cry__getAtom2Num();
atom2num_len = size(atom2num,1);

% remove white space
atomicsymbol = strtrim(atomicsymbol);

% search for the cell array index that gives us the atomic symbol 
found = false;
for index = 1:atom2num_len
    if(strcmp(atom2num(index, 1), atomicsymbol))
        found = true;
        atomicnumber = cell2mat(atom2num(index, 2));
        break;
    end
    
end
if found == false
    error(sprintf('Atomic number for atomic symbol %s not found.\nCheck routine cry__getAtom2Num().', atomicsymbol))
end


