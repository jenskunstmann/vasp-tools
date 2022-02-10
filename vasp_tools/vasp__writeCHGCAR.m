function vasp__writeCHGCAR(filename, structure, dens)
% generate a VASP CHGCAR file
%
% we save the crystal structure and density in two separate files and cat
% them together afterwards
% we save the charge part with one call of 'save' because calling fprintf()
% npnts times takes forever (5 min for a 150 MB file), while save needs
% only 12 sec for the same data
%
% we want to write out 5 points per line, as some codes rely on this
% =1, impies that the output format is exactly as in VASP, SLOW! 
% DO NOT USE THIS OPTION: PART OF THE CODE IS BUGGY
%matchformat = 0;

% write POSCAR file 'cryst'
vasp__writeCONTCAR('cryst', structure)

% dimension of the charge density mesh
[nx,ny,nz] = size(dens);
npnts = nx*ny*nz;

densvec = reshape(dens, npnts, 1);   % transform into a vector

% THIS PART OF THE CODE IS BUGGY, AS IT MESSES UP THE DENSITY COMPLETELY
% transform the density into a (nlines X pointsperline) matrix; add zeros
% to the last line, if such a matrix does not trivially exist
% pointsperline = 5;         % 5 points per lines are the VASP default
% nlines = floor(npnts/pointsperline);  % number of lines with 'pointsperline' data poins
% nmissing = npnts - pointsperline*nlines;  % number "missing" data points if only nlines are written
% if nmissing > 0    % fill up vector swith zeros to enable reshape() to work    
%     nadd = pointsperline - nmissing;
%     densvec(npnts+1:npnts+nadd) = zeros(nadd,1); % add a few data points
%     nlines = nlines + 1;    
% end
% densvec = reshape(densvec, nlines, pointsperline);   % transform into a matrix

% add newline and size of the data grid to the crystal structure
fh = fopen('cryst', 'a');
fprintf(fh, '\n %4d %4d %4d\n', nx,ny,nz);

% now ...
display(sprintf('Writing charge density'))
% if matchformat
%     % write the field line by line
%     for line = 1:nlines
%         fprintf(fh, ' %4.11E', densvec(line,:));
%         fprintf(fh, '\n');
%     end
%     fclose(fh);
%     movefile('cryst', filename)
    
%else
    % close the 'cryst' file
    fclose(fh);
    
    % this is a fast and efficient method to save the big field
    save('charge', '-ascii', '-double', 'densvec');     

    % merge crystal structure and charge density file
    commandstr = sprintf('cat > %s cryst charge', filename);
    [status, result] = system(commandstr);

    % remove the temporary files: cryst charge
    delete cryst charge
%end