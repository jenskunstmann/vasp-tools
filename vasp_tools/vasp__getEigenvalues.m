function [kpnt_pos, eval_up, eval_down] = vasp__getEigenvalues(varargin)
% use global vaiable SYS to read in eigenvalues either from EIGENVAL or
% OUTCAR file, in case of scf band structure calculations, read in IBZKPNT
% file and remove block of scf k-points 
%
% kpnt_pos(kpoint, component): position of the k-point in reduced
%                              coordinates,  i.e. in units of the 
%                              reciprocal lattice vectors
% eval_up/down(kpoint, band) : band energies at a particular k-point 

global SYS

% define file source by arguments
switch(nargin)
    case(0)
        % make filename
        file_eigenval = sprintf('%s/%s',SYS.path,SYS.eigenval)
        
    case(1)
        % take filename
        file_eigenval = varargin{1}
        
    otherwise
        error('Wrong number of function arguments.')
end



% figure out if eigenvalues are read from EIGENVAL or OUTCAR
indices = strfind(file_eigenval, 'EIGENVAL');
if isempty(indices)
    % read from OUTCAR
    readeigenval = false;
else
    % read from EIGENVAL
    readeigenval = true;
end

% read EIGENVAL/OUTCAR file, if required remove the block of scf k-points 
if (isempty(SYS.ibzkpt))    % traditional NON-scf band calc    
    if readeigenval
        [kpnt_pos, eval_up, eval_down] = vasp__readEIGENVAL(file_eigenval); 
    else
        [kpnt_pos, eval_up, eval_down] = vasp__readOUTCAR(file_eigenval, 0);
    end
    
else                        % scf band calc
    file_ibzkpnt  = sprintf('%s/%s',SYS.path,SYS.ibzkpt)
    nkpnts = vasp__readKPOINTS(file_ibzkpnt);

    if readeigenval
        [kpnt_pos, eval_up, eval_down] = vasp__readEIGENVAL(file_eigenval, nkpnts); 
    else
        [kpnt_pos, eval_up, eval_down] = vasp__readOUTCAR(file_eigenval, nkpnts);
    end        
end