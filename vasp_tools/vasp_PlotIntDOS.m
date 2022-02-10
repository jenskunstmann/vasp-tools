function vasp_PlotIntDOS()
% plot integrated DOS from VASP
% works for spin-polarized calculations

global SYS

% read DOSCAR file
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar);
[e ddos idos] = vasp__readDOSCAR(file_doscar);

% define the 'dos' structure
dos.e = e;
dos.ddos = [idos]; % 0.9*idos];
dos.lstyle = {'-','-'};     % continuous lines
dos.legend = {};       % no legends

% set dos.sort
% non-spin-polarized if 'ddos' has only one column
if (size(dos.ddos,2) == 1)
    % non spin-polarized, single-panel plot
    dos.sort   = {[1]};  
    dos.color = [0 0 0]; 
else
    % spin-polarized
    % plot (spin-up) and (spin-down) curve in one same panel      
    dos.sort   = {[1 2]};    
    dos.color = [1 0 0; 0 0 1]; % spin-up = red, spin-down=blue
end

dos.plotinset = 1;        % true:1: plot inset, false:0: plot only DOS

% plotting ranges
dos.emin = min(e); 
dos.emax = max(e); 
dos.dosmax = max(idos(:,1)); 
dos.dosmin = 0;    

% plotting ranges of the INSET
dos.inset_emin = -3;
dos.inset_emax = 3; 
dos.inset_dosmax = dos.dosmax;
dos.inset_dosmin = 0; 

% shift inset and main plot in y direction by insetshift*dosmax
dos.insetshift = 1.0;
%%% END: user settings %%%

[mainax, insetax] = vasp__plotDOS(dos);
title(mainax, SYS.ID,'Interpreter','none');
ylabel(mainax, 'iDOS (states / atom)');








