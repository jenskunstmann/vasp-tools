function vasp_NNeigborAnalysis()
% show the nearest neigbors and bond lengths of an atomic structure
% results depend on the specified cutoff

global SYS

%%%% radius of atomic sphere to search neigbors, within
cutoff = 4;
%%%%%%

file_contcar = sprintf('%s/%s', SYS.path, SYS.contcar)
%file_contcar = '/exports-hal/users/home/jkunstmann/res/GaAs/structure/POSCAR';
crystal = cry__readCONTCAR(file_contcar);

cry__NNeigborAnalysis(crystal, cutoff)

% % nearest neigbors
% nntable = cry__FindNearestNeighbors(crystal, cutoff);
% for at = 1:natoms
%     display(['neigbors of atom ' num2str(at)])
%     nntable{at}.atomID
%     clear nndist
%     for nn = 1:size(nntable{at}.distvec,1)
%         nndist(nn) = norm(nntable{at}.distvec(nn,:));    
%     end
%     nndist
%     %nntable{at}.distvec
% end
