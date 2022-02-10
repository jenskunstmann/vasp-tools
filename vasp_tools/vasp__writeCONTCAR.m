function vasp__writeCONTCAR(filename, structure)
% write the CONTCAR file
% 
% structure.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% structure.types(type)       = number of atoms per type 'type'
% structure.atompos(at,comp)  = atomic coordinates as given in the file, i.e. either
%                     in Cartesian or in direct coordinates
% structure.comment      = first line of the file (comment)
% structure.format       = format line, i.e. 'Direct', 'Cartesian', etc. 
% structure.vasp52       = true/false; which file format
% structure.atsymb       = line with the atomic symbols if vasp52=true

%%%%% write out data in CONTCAR file format
% open file
display(sprintf('Writing crystal structure to file %s', filename))
fh = fopen(filename, 'w');

fprintf(fh, '%s\n   1.00000000000000\n', structure.comment);
% general format
% fprintf(fh, '   %4.16f %4.16f %4.16f\n', structure.latt(1,:));
% fprintf(fh, '   %4.16f %4.16f %4.16f\n', structure.latt(2,:));
% fprintf(fh, '   %4.16f %4.16f %4.16f\n', structure.latt(3,:));

% this is the CHGCAR VASP format, required exactly by the 'Bader' code to
% work
fprintf(fh, '  %11.6f %11.6f %11.6f\n', structure.latt(1,:));
fprintf(fh, '  %11.6f %11.6f %11.6f\n', structure.latt(2,:));
fprintf(fh, '  %11.6f %11.6f %11.6f\n', structure.latt(3,:));


% write the atom type line(s)
if structure.vasp52 == true
    fprintf(fh, ' %s\n', structure.atsymb);
end
for i=1:size(structure.types,1)
    fprintf(fh, ' %3d', structure.types(i));
end

% format line
fprintf(fh, '\n%s\n', structure.format);

% atomic positions
for i=1:size(structure.atompos,1)
    %fprintf(fh, '   %4.16f %4.16f %4.16f\n', structure.atompos(i,:));
    fprintf(fh, ' %9.6f %9.6f %9.6f\n', structure.atompos(i,:));
end

% close file
fclose(fh);

