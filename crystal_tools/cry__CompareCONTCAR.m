function cry__CompareCONTCAR(filename1, filename2)

% read in files
crystal1 = cry__readCONTCAR(filename1);
crystal2 = cry__readCONTCAR(filename2);

% compare
cry__CompareStructures(crystal1, crystal2)