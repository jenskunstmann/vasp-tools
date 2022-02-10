function nnmax = cry__GetMaximumNumberOfNearestNeighbors(nntable)
% determines the maximum number of NN that is present

nnmax = 0;
natoms = length(nntable);       % number of atoms

for at = 1:natoms
    numNN = length(nntable{at}.atomnum); % number of nn of this atom
    if numNN > nnmax
        nnmax = numNN;
    end
end
