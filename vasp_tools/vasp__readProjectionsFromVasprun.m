function  vasp__readProjectionsFromVasprun()
% generate PROCAR.mat file from vasprun.xml
% in the xml file the band characters are given with 5 significatn digits
% while in PROCAR there are only 4. This can be a problem in big systems.
%
% Because vasprun.xml for large systems is a very huge file, it cannot be
% parsed easily with standard xml routines (not enough memory). It is
% therefore read in line by line.
%
% kpnt_pos(kpoint, component): 
%         position of the k-point in reduced coordinates,  i.e. in units of
%         the reciprocal lattice vectors 
% eval(kpoint, band) : 
%         band energies at a particular k-point 
% bandchar(kpnt,band,atom,orbital): 
%         the orbital=(lm) and atom projected character of a band

%%%% user data
lorbit = 11;        % = 10 (l projections), = 11 (lm projections)
xmlfile = '/data/jk/tmp/bands-4/vasprun.xml'
matfile = '/data/jk/tmp/bands-4/PROCAR.mat';
%%%%

% open vasprun.xml file
fh = fopen(xmlfile, 'r');

% tell user 
display('Reading projections from vasprun.xml file -> PROCAR.mat')

% read in k-points
display('reading k-points')
FindInFile(fh, '<kpoints>');      
FindInFile(fh, '<varray name="kpointlist" >');     
% do it
nkpnts = 0;
while(true)
    [Q,l] = fscanf(fh, ' <v>  %g %g %g </v> ', 3);
    if l == 3   % line sucessfully read
        nkpnts = nkpnts + 1;
        kpnt_pos(nkpnts,:) = Q';
    else
        break;
    end                     
end
nkpnts

FindInFile(fh, '<parameters>');      
%FindInFile(fh, '<separator name="electronic" >');      
tline = FindInFile(fh, 'NBANDS');
[Q,l] = sscanf(tline, ' <i type="int" name="NBANDS"> %i </i> ', 1);
if l ~= 1
     error('error reading file');   
end      
nbands = Q    

FindInFile(fh, '<atominfo>');         
[Q,l] = fscanf(fh, ' <atoms> %i </atoms> ', 1);
if l ~= 1
     error('error reading file');   
end      
natoms = Q                      

% find projections
display('reading eigenvalues')
FindInFile(fh, '<projected>');       
% skip header lines
skipline(fh, 9);   

% read eigenvalues; first child of <projected>
eval = zeros(nkpnts, nbands);
for kpnt = 1:nkpnts
    skipline(fh, 1); % <set comment="kpoint #">

    % read band energy
    [Q,l] = fscanf(fh, ' <r>  %g %g </r> ', [2 nbands]);
    if l ~= 2*nbands
         error('error reading file');   
    end            
    eval(kpnt,:) = Q(1,:);   

    skipline(fh, 1);  % </set>
end

% read in band characters;
% next child of <projected>
display('reading projections')
FindInFile(fh, '<array>');
switch(lorbit)
    case 10
        skipline(fh, 9);
    case 11
        skipline(fh, 15);
end

bandchar = zeros(nkpnts, nbands, natoms, 9);
for kpnt = 1:nkpnts
    display(['  reading k-point ' num2str(kpnt)])
    skipline(fh, 1); % <set comment="kpoint #">

    for band = 1:nbands
        skipline(fh, 1); % <set comment="band #">
        
        switch(lorbit)
            case 10
                % read characters: l projections 
                [Q,l] = fscanf(fh, ' <r> %g %g %g </r> ', [3 natoms]);
                if l ~= 3*natoms
                     error('error reading file');   
                end            
                bandchar(kpnt,band,:,1) = Q(1,:)'; % s 
                bandchar(kpnt,band,:,2) = Q(2,:)'; % p
                bandchar(kpnt,band,:,5) = Q(3,:)'; % d                              
                
            case 11
                % read characters: lm projections 
                [Q,l] = fscanf(fh, ' <r> %g %g %g %g %g %g %g %g %g </r> ', [9 natoms]);
                if l ~= 9*natoms
                     error('error reading file');   
                end            
                bandchar(kpnt,band,:,:) = Q';

        end
        skipline(fh, 1);    % </set>    
    end

    skipline(fh, 1); % </set>
end               

% close file
fclose(fh);

% save data in mat-file, v7.3 is required for very big arrays
display(['writing file: ', matfile])
save(matfile,'bandchar','eval','kpnt_pos', '-v7.3');
    


function tline = FindInFile(fh, pattern_str)
% read file line by line and find next occurance of 'pattern_str'
% current position in file is contained in 'fh'

found = false;

tline = fgetl(fh);
while ischar(tline)        
    if strfind(tline, pattern_str)
        found = true;
        break
    end
    tline = fgetl(fh);
end

if ~found    
    error('search pattern not found in file')
end


          