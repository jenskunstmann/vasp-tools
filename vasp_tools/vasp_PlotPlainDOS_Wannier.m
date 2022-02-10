function vasp_PlotPlainDOS_Wannier(varargin)
% plot DOS from VASP
% works for spin-polarized calculations

global SYS

% allow to chage the Fermi energy
switch(nargin)
    case(1)
        DeltaEF = varargin{1};      
        
    otherwise
        DeltaEF = 0;
end

% read DOSCAR file
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar)
[e ddos idos natoms] = vasp__readDOSCAR(file_doscar); 
efermi =  vasp__getEFermi(file_doscar) 
efermi_shifted = efermi + DeltaEF;

% define the 'dos' structure
eshift = 0;
dos.e = e - eshift;
dos.ddos = ddos;
dos.lstyle = {'-','-'};     % continuous lines
dos.legend = {};       % no legends

% dos.ddos = [ddos ddos 2*ddos 0.5*ddos];
% dos.lstyle = {'-','-','-','-'};     % continuous lines
% dos.legend = {'DOS','DOS','2DOS','0.5DOS'};       % no legends

% set dos.sort
% non-spin-polarized if 'dos' has only one column
if (size(dos.ddos,2) == 1)
    % non spin-polarized, single-panel plot
    dos.sort   = {[1]};  
    dos.color = [0 0 0]; 
    
    % show the band gaps
    display('  E_start    E_end    Delta E_gap')
    [estart, eend] = el__getBandGaps(e, ddos);
    gaps = [estart'  eend' eend' - estart']
else
    % spin-polarized
    % plot first curve of ddos (spin-up) in the 
    % uppper and the second curve (spin-down) in the lower panel     
    dos.sort   = {[1],[2]};    
    dos.color = [1 0 0; 0 0 1]; % spin-up = red, spin-down=blue
%     dos.sort   = {[1 3],[2 4]};    
%     dos.color = [1 0 0; 0 0 1; 0 1 0; 0 1 1]; % spin-up = red, spin-down=blue

    % show the band gaps
    display('  E_start    E_end    Delta E_gap')
    [estart, eend] = el__getBandGaps(e, ddos(:,1));
    gaps1 = [estart'  eend' eend' - estart']
    
    [estart, eend] = el__getBandGaps(e, ddos(:,2));
    gaps2 = [estart'  eend' eend' - estart']
end

dos.plotinset = 0;        % true:1: plot inset, false:0: plot only DOS

% plotting ranges
dos.emin = min(e); % SYS.emin; %
dos.emax = max(e); % SYS.emax; % 
dos.dosmax = max(ddos(:,1));
dos.dosmin = -dos.dosmax;    % only significant for spin-pol. calcs

% plotting ranges of the INSET
dos.inset_emin = -3;
dos.inset_emax = 3; 
dos.inset_dosmax = 0.01*dos.dosmax;
dos.inset_dosmin = dos.dosmin; % only significant for spin-pol. calcs

% shift inset and main plot in y direction by insetshift*dosmax
dos.insetshift = 1.0;
%%% END: user settings %%%

[mainax, insetax] = vasp__plotDOS(dos);

hold on
%%%%%%%%%
%read in the DOS vom wannier90 and plot it
%natoms = 4;  % <- get is automatically
file_dos_wann = sprintf('%s/%s/%s', SYS.path, SYS.wandir, 'wannier90-dos.dat')
dat_wann = load(file_dos_wann);
e_wann  = dat_wann(:,1) - efermi_shifted;   % x-axis = energy
dos_wann = dat_wann(:,2)/natoms;            % DOS is normalized to number of atoms
plot(e_wann, dos_wann, 'r')
[estart, eend] = el__getBandGaps(e_wann, dos_wann);
gaps = [estart'  eend' eend' - estart']


title(mainax, SYS.ID,'Interpreter','none');
ylabel(mainax,'DOS (states / eV atom)');
xlabel(mainax, 'E - E_F (eV)');
SetFontsInFigure(20);







