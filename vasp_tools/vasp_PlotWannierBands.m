function vasp_PlotWannierBands(varargin)
% make a plain bandstructure plot for wannier90 output, 
% k-points and distances correspond to Carestian coordinates in reciprocal space
% our own k-labels are used
%
% USAGE: vasp_PlotWannierBands([DeltaEF])
%
% DeltaEF = change of Fermi energy, allows to align different band plots

global SYS

% allow to chage the Fermi energy
switch(nargin)
    case(1)
        DeltaEF = varargin{1};      
        
    otherwise
        DeltaEF = 0;
end

emin = SYS.emin;        % plotting range
emax = SYS.emax;                
klabels = SYS.klabels;          % special point labels
linespec = 'k';  

% get Fermi energy from DOSCAR
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar)
efermi =  vasp__getEFermi(file_doscar) 
efermi_shifted = efermi + DeltaEF;
%efermi = 0


% wannier90 files
file_wout = sprintf('%s/%s/%s', SYS.path, SYS.wandir, 'wannier90.wout')
file_eval = sprintf('%s/%s/%s', SYS.path, SYS.wandir, 'wannier90_band.dat')
file_gnu = sprintf('%s/%s/%s', SYS.path, SYS.wandir, 'wannier90_band.gnu')

% read in these files
nwannier = vasp__readWannier90Wout(file_wout);  % # wannier functions = # bands
gnu_str = fileread(file_gnu); 
%contstrlen = length(gnu_str);
dat_wann = load(file_eval);
[len,~] = size(dat_wann);               % len X 2 matrix, with (hidden) nwannier blocks

% plot the bands
e_wann = dat_wann(1:(len/nwannier),1);    % x-axis = energy
eval_wann = reshape(dat_wann(:,2), [len/nwannier nwannier])-efermi_shifted; % eigenvalues 
plot(e_wann, eval_wann, linespec)

hold on 

% find the positions of the k-labels in the gnu file
indx_start = regexp(gnu_str, '(');
indx_end  = regexp(gnu_str, ')');
line = gnu_str(indx_start:indx_end);  % string between '( )'
[start, ende] = regexp(line, '[-+]?[0-9]*\.?[0-9]+'); % position of numbers
nlabels = length(start);  % number of numbers/labels
for i = 1:nlabels
    label_pos(i) = str2num(line(start(i):ende(i)));
end
xmax = label_pos(nlabels);


% plot labels and lines
% plot zeroline = EF
plot([0 xmax], [0 0], '--k');

% draw vertical separation lines
vtmp = ones(1,nlabels);                             % temporary row vector
xpnts = [label_pos; label_pos];
ypnts = [emin*vtmp; emax*vtmp];         % 2x(nlabels) dimensional matrix of y positions
plot(xpnts, ypnts, '-k');  

% plot special point labels
elabel = emin - (emax-emin)*.05;    % y position of the labels
text(label_pos, elabel*vtmp, klabels, ...
    'HorizontalAlignment', 'center', 'Interpreter','latex');

% format the plot
set(gca, 'XTick', []);      % remove x axis ticks
ylabel('E-E_F (eV)');      % y label
box on                      % display boxed axes
axis([0 xmax emin emax]);   % visible range
title(SYS.ID,'Interpreter','none');
SetFontsInFigure(20);

% save band energies in an ASCII file, can be read in xmgrace as block data
%ematrix = [x' eval_up - efermi];
%save('bands.dat', 'ematrix', '-ASCII')