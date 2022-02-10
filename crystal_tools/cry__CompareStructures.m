function cry__CompareStructures(crystal1, crystal2)
% compares two crystal structues and calculates the mean square
% displacement per Cartesian component

natoms = length(crystal1.atomnum);

if natoms ~= length(crystal2.atomnum)
    error('The two crystals contain a different number of atoms.')
end

% calculate the differnce of the positions
diff = crystal1.atompos - crystal2.atompos
sumerrsquare = sum(diff .* diff, 1);
MeanSquareDisplacement = sqrt(sumerrsquare/natoms)

