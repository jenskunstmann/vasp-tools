function [kpnt_pos, eval_up, eval_down] = vasp__readOUTCAR(filename, nkpntstoskip)
% read eigenvalues from OUTCAR instead of EIGENVAL
%
% syntax: [kpnt_pos, eval_up, eval_down] = vasp__readOUTCAR(filename, nkpntstoskip)
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

% NO SPIN POLARIZATION IMPLEMENTED YET!!!
eval_down = [];

% read file content into string
contstr = fileread(filename);
contstrlen = length(contstr);

% get some general parameters
efermi = getvalue(contstr, 'E-fermi :');
nkpnts = getvalue(contstr, 'NKPTS =');
nbands = getvalue(contstr, 'NBANDS=');


%%%% read k-point coordinates 
searchstr = ' k-points in reciprocal lattice and weights: k-points along high symmetry lines';
indices = strfind(contstr, searchstr);
startindex = indices(1) + length(searchstr);  

[A, count] = sscanf(contstr(startindex:contstrlen), ' %g %g %g %g ', [4 nkpnts]);
if count ~= 4*nkpnts
    error('error reading file');         
end;     

kpnt_pos = A(1:3, (nkpntstoskip+1):nkpnts)';
   
 
 
%%%% read eigenvalues
searchstr = 'band No.  band energies     occupation';
indices = strfind(contstr, searchstr);
if nkpnts ~= length(indices)
    error('OUTCAR does not contain all k-points.')
end

for kpnt = 1:nkpnts
    startindex = indices(kpnt) + length(searchstr); 
     
    [A, count] = sscanf(contstr(startindex:contstrlen), ' %g %g %g ', [3 nbands]);
    if count ~= 3*nbands
        error('error reading file');         
    end;     
    
    eval_up_tmp(kpnt,:) = A(2,1:nbands);    
end

% skip first block of k-points if nkpntstoskip > 0
eval_up  = eval_up_tmp((nkpntstoskip+1):nkpnts,:);



function value = getvalue(str, searchstr)
% return numerical value after keyword

index = strfind(str, searchstr);
if isempty(index)
    error(sprintf('keyword <%s> not found.', searchstr))
end

startindex = index(1) + length(searchstr);
endindex   = startindex + 50;
[value, count] = sscanf(str(startindex:endindex), ' %f ');
if count ~= 1
    error('error reading file');         
end;