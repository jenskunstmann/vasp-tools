function vasp_PlotUnfoldedBands()
% read in the band structure file produced by the unfolding tool BandUP and
% plot the effective band strucure (ebs) produced by unfolding
%
% run vasp_Start() first and load a system where unfolding was done before
%
% the end point of one k-line and the start-point of the next one are the
% same, thus the data point where two k-lines meet is redundant, however I
% checked that the surf() function does not care and the spectral function
% is the same no matter whether we keep or remove the redundant point
%
% in BandUP the spectral function A(k,E) at any k-point is integrated over
% the energy range E + dE  
% dN(k,E) = Int_(E)^(E+dE) A(k,E') dE' 
% = the number of primitve cell (PC) bands (states), at the PC wave vector
% k in that energy interval    
%
% dN(k,E) is what is actually calculated and plotted
%
% HINT:
% dE >= 0.2 for good plotting
% dE = 0.001 for good band mass fitting


% input parameters: k-labels and filename
global SYS
klabels = SYS.klabels;    % e.g. {'M', '$\Gamma$', 'K', 'M'};
file_ebs = sprintf('%s/%s',SYS.path,SYS.ebs)

% open file
fh = fopen(file_ebs, 'r');

% skip header lines and determine the number of columns in the data part
% there are two formats: with and without spin giving 3 or 5 columns, resp.
% however the first 3 columns are always the same
ncolumns = skipheader(fh);

% read in all data points
[Q,l] = fscanf(fh, ' %g ');
if l < ncolumns
    error('error reading file');                 
end;   

fclose(fh);

% determine the dimensions of our matrices
npnts = l/ncolumns;                         % number of data points
data = reshape(Q, ncolumns, npnts)';        % data as in the file as npnts x ncolumns matrix
npnts_e = npnts - nnz(data(:,1));    % number of energy points = number of zero-valued elements in the x-vector
npnts_k = npnts/npnts_e;             % number of k-points 

% generate the 2D matrices X,Y and intensity I for plotting
X = reshape(data(:,1), npnts_e, npnts_k);
Y = reshape(data(:,2), npnts_e, npnts_k);
I = reshape(data(:,3), npnts_e, npnts_k);

x = X(1,:);   % x-axis (k-axis)
y = Y(:,1)';    % y-axis (E-axis)


% min/max values for k (x), E (y) / plotting range
limk = [min(x) max(x)];
%limE = [SYS.emin SYS.emax]; 
limE = [min(y) max(y)];

% plot the effective band structure spectral function dN(k,E)
hold on
obj = surf(X, Y, I, ...
         'EdgeColor','none','BackFaceLighting','lit',...
         'AmbientStrength',0.5, 'FaceAlpha', 1);

% add/remove interpolation of the data points
set(obj, 'FaceColor','interp')         
         
% make colobar
% read in a previously generated a color map
% mycmap = get(gcf, 'Colormap');
colorbar 
load('MyColormaps', 'mycmap')
set(gcf, 'Colormap', mycmap)

% plot zeroline = EF
plot(limk, [0 0], '--k');


% draw vertical separation lines
% the NUMBER of k-labels is used to determine the position of the
% separation lines because VASP always uses the SAME number of k-points
% along each k-line
nlabels = size(klabels,2);
xind = [1 npnts_k/(nlabels-1)*(1:nlabels-1)];   % indices 'xind' of  x(xind) that correspond to the posn. of a k-point label
vtmp = ones(1,nlabels);                             % temporary row vector
xpnts = [x(xind); x(xind)];             % 2x(nlabels) dimensional matrix of x positions
ypnts = [limE(1)*vtmp; limE(2)*vtmp];         % 2x(nlabels) dimensional matrix of y positions
plot(xpnts, ypnts, '-k');  

% plot special point labels
elabel = limE(1) - (limE(2)-limE(1))*.05;    % y position of the labels
text(x(xind), elabel*vtmp, klabels, ...
    'HorizontalAlignment', 'center', 'Interpreter','latex');

% axes settings
set(gca, 'XTick', []);      % remove x axis ticks
ylabel('E-E_F (eV)');      % y label
xlim(limk)
ylim(limE)
view(0, 90)        % camera position
shg;                % bring to front
box on             % axes boxed

title(SYS.ID, 'Interpreter', 'none');
SetFontsInFigure(20);


  
function ncolumns = skipheader(fh)
% skip header lines starting with '#' signs
% and further determine the number of columns in the data part

% first find out how many lines of header there are
hlines = 0;     % number of header lines
while true
    head = fgetl(fh);
    %disp(head)
    
    if head(1) ~= '#'
        break;
    else
        hlines = hlines + 1;
    end        
end;

[~, ncolumns] = sscanf(head, ' %g ');

% problem is that the first line of data was read in and is gone now;
fseek(fh, 0, 'bof' );    % return to the beginning of the file
skipline(fh, hlines);   % and skip only 'hlines' 

