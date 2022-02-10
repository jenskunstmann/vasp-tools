function [wkpnt_pos, weval_up, weval_down] = vasp__readEIGENVAL(filename, varargin)
% read the band structure data from an EIGENVAL file of VASP
% for non-spin-polarized calcs. 'eval_down=0'
%
% syntax: [kpnt_pos, eval_up, eval_down] = vasp__readEIGENVAL(filename [,nkpntstoskip])
% 
% for hybrid functionals, GW, etc. the bands calculation is fully
% self-consistent. The k-lines follow after a block of regular equally
% spaced k-points. then the first 'nkpntstoskip' k-points, that do not
% correspond to k-lines, can be removed.
%
% kpnt_pos(kpoint, component): position of the k-point in reduced
%                              coordinates,  i.e. in units of the 
%                              reciprocal lattice vectors
% eval_up/down(kpoint, band) : band energies at a particular k-point 

% PROBLEMS: the EIGENVAL format changed to include a column with the
% occupation numbers. Currently there is no backward compatible implementation


% define file source by arguments
switch(nargin)
    case(2)
        nkpntstoskip = varargin{1};      
        
    otherwise
        nkpntstoskip = 0;
end

% open file
fh = fopen(filename, 'r');

% read first line
[Q,l] = fscanf(fh, ' %i  %i  %i %i ', [4 1]);
if l ~= 4
    error('error reading file: first line');         
end;

% if the 4rth number is '2' then the calc was spin-polarized
if Q(4) > 1
    spinpol = 1;
else
    spinpol = 0;
    eval_down = [];
end;

% skip 4 lines of header
skipline(fh, 4);

% read dimension line 
[Q,l] = fscanf(fh, ' %i  %i  %i ', [3 1]);
if l ~= 3
    error('error reading file, dimension line');         
end;

nkpnts = Q(2);   % number of k-points
nbands = Q(3);   % number of bands (eigenvalues) per k-point

% determine version of file
skipline(fh, 1);
if(spinpol)
    error('Not working yet. Please implement the version check for spinpolarized calculations!')
else
    % EIGENVAL file changed to contain 3 columns, older versions
    % contained 2 columns;
    % VASP uses fixed format, so if the band number > 999 it turns into
    % '***', to make this work we consider it as a fixed column string
    % and do not read it in; it must be fixed column format to work in
    % all cases    
    line = fgetl(fh);    
    [Q,ncolumns] = sscanf(line, ' %*3s %g %g ', 2 );       
    if ncolumns == 1 
        formatstr = ' %*3s %g ';    % older format        
    elseif ncolumns == 2
        formatstr = ' %*3s %g %g '; % recent format        
    else
        error('reading EIGENVAL: unknown number of columns.')        
    end
end
% return to the beginning of the file 
frewind(fh)
skipline(fh, 7);

% loop over all kpoints
for kpnt = 1:nkpnts
    
    % read position of k-point
    [Q,l] = fscanf(fh, ' %g %g %g %g ', [4 1]);
    if l ~= 4
        error('error reading file, kpoints');                 
    end;   
    kpnt_pos(kpnt,:) = [Q(1) Q(2) Q(3)];
    
    % read block of 'nbands' eigenvalues
    if(spinpol)
        % read two energies if spin-polarized
        % !!!! not updated for newer versions   !!!
        [Q,l] = fscanf(fh, ' %*3s %g %g ', [2 nbands]);
        if l ~= 2*nbands
            error('error reading file, band energies');                 
         end;
             
        eval_up(kpnt,:) = Q(1,1:nbands);
        eval_down(kpnt,:) = Q(2,1:nbands);         
        
    else
        % non-spin-polarized
        [Q,l] = fscanf(fh, formatstr, [ncolumns nbands]);     
        if l ~= ncolumns*nbands
            error('error reading file, band energies');                 
        end         
        eval_up(kpnt,:) = Q(1,:); 
        
    end                    
end

% close file
fclose(fh);

% skip first block of k-points if nkpntstoskip > 0
wkpnt_pos = kpnt_pos((nkpntstoskip+1):nkpnts,:);
weval_up  = eval_up((nkpntstoskip+1):nkpnts,:);

if(isempty(eval_down))     % non-spin-polarized if 'eval_down' is empty
    weval_down = [];
else    
    weval_down = eval_down((nkpntstoskip+1):nkpnts,:);
end


return;


