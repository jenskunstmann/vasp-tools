function kpoints = el__MakeKLines(specpoints, subdiv)
% generate the kpoints of the K-lines for band structures
% boundary points are generated twice (to be compatible with VASP);
% different the total number of points = (numspecpoints-1)*(subdiv+1)
% 
% specpoints(:,comp) = coordinates of the special points in units of the rlv
% subdiv             = number of subdivisions along one k-line
%
% kpoints(:,comp)    = coordinates of the k-points in units of the rlv
%
% the number of points per line = subdiv + 1; start and endpoints of
% subsequent k-lines coincide

numlines = size(specpoints,1) - 1;

x = [0:(1/subdiv):1]';
kpoints = zeros((numlines*(subdiv+1)), 3);
maxindx = 0;
for line = 1:numlines
    pstart = specpoints(line,:);
    pend   = specpoints(line+1,:);    
    dir = pend - pstart;
    
    minindx = maxindx + 1;          % no +1 as atart and end points of different lines are the same
    maxindx = minindx + subdiv;    
    kpoints(minindx:maxindx,:) = pstart(ones(1,subdiv+1),:) + x*dir; 
    
end
