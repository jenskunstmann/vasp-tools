function cry__NNeigborAnalysis(crystal, cutoff)
% show the nearest neigbors and bond lengths of an atomic structure
% results depend on the specified cutoff
%
% crystal = 'crystal' data structure
% cutoff = radius of atomic sphere to search neigbors, within

natoms = length(crystal.atomnum);

% nearest neigbors
nntable = cry__FindNearestNeighbors(crystal, cutoff);
for at = 1:natoms
    display(['neigbors of atom ' num2str(at)])
    nntable{at}.atomID
    clear nndist
    for nn = 1:size(nntable{at}.distvec,1)
        nndist(nn) = norm(nntable{at}.distvec(nn,:));    
    end
    nndist
    %nntable{at}.distvec
end
