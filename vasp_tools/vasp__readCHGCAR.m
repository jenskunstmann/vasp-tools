function [crystal, density] = vasp__readCHGCAR(filename)
% read in VASP CHG/CHGCAR/ELFCAR/.. files
%
% filename = name of CHGCAR file  
%
% atomic structure defined by
% crystal.origin(comp) = Cartesian coordinates of the origin of the
%                        structure
% crystal.latt(numvec,comp) = Bravais lattice vectors; numvec=1,2,3; comp=1,2,3
%                     Cartesian component
% crystal.atompos(atomID,comp)  = Cartesian coordinates of the atom; comp=1,2,3
% crystal.typeatomnum(atomID)       = atomic number of the atom 
%
% density(nx,ny,nz)    = density mesh

%
% current limitations: 
% see cry__readCONTCAR()
% MIND: the field values are NOT rescaled. charge density needs to be
% divided by 1/V_cell

% open file
fh = fopen(filename, 'r');

% read crystal structure
%[structure] = vasp__readCONTCAR(fh);
[crystal] = cry__readCONTCAR(fh);

% read the size of the data grid
[Q,l] = fscanf(fh, ' %f  %f  %f  ', [3 1]);
if l ~= 3
    error('error reading file');         
end;
nx = Q(1);
ny = Q(2);
nz = Q(3);

% read data block
% npnts = number of data points
npnts = nx*ny*nz;
[densvec,l] = fscanf(fh, '%f', [npnts 1]);
if l ~= npnts
    error('error reading file');         
end;

% reshape the density vector into a matrix of the actual dimensions
density = reshape(densvec, nx, ny, nz);

% close file
fclose(fh);


% OLD CODE TAKEN FROM ABOVE
% % read comment line 
% more.comment = fgetl(fh);
% 
% % read the scaling factor
% [scale,l] = fscanf(fh, ' %f ', 1);
% if l ~= 1
%     error('error reading file');         
% end;
% 
% % read the lattice vectors
% [latt,l] = fscanf(fh, ' %f  %f  %f  ', [3 3]);
% if l ~= 9
%     error('error reading file');         
% end;
% latt = latt';
% 
% % rescale lattice vectors
% if scale > 0
%     latt = latt * scale;
% else
%     error('lattice volume scaling not yet implemented.')
% end
% 
% % read atom type line(s)
% line = fgetl(fh);
% [types,ntypes] = sscanf(line, ' %d %d %d %d %d %d %d %d %d %d ', 10);
% if ntypes < 1           
%     % if reading was not successfull we have VASP 5.2 file format
%     more.vasp52 = true;    
%     more.atsymb = line;  % line with atomic symbols
%     % read next line = atom type line
%     line = fgetl(fh);
%     [types,ntypes] = sscanf(line, ' %d %d %d %d %d %d %d %d %d %d ', 10);
%     if ntypes < 1 
%        error('error reading file');         
%     end    
% else
%     % if reading was successfull we have the old VASP file format
%     more.vasp52 = false;
% end;
% natoms = sum(types);
% 
% % read format line; 'selective dynamics' will not work, yet
% more.format = fgetl(fh);
% 
% % read the atomic positions; 
% % reading them in as a (natoms X 3) matrix does not work for some reason
% % therefore we read it in as (3 X natoms) matrix and transpose afterwards
% [atompos,l] = fscanf(fh, ' %f  %f  %f ', [3 natoms]);
% if l ~= natoms*3
%     error('error reading file');         
% end;
% atompos = atompos';
% 
% % conversion from the old format, done in order not to change the code too
% % much
% structure.latt     = latt;
% structure.types    = types;
% structure.atompos  = atompos;
% structure.comment  = more.comment;
% structure.format   = more.format;
% structure.vasp52   = more.vasp52;
% structure.atsymb   = more.atsymb;



