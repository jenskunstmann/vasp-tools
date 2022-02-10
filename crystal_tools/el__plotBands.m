function [x, ax] = el__plotBands(bands)
% make a (plain/projected) bandstructure plot as specified by the 'bands'
% structure
%
% bands.emax = plotting range
% bands.emin 
% bands.klabels{} =  special point labels
% bands.eval(kpnt,band) = band energies with Fermi level at E=0
% bands.kpnt_pos(kpnt,:) = positions of kpoints (in any given unit)
% bands.linespec = 'linespec' of the lines of the bands
% bands.color(:) = color of the lines
% bands.bchar(kpnt,band) = band characters = fatbands
% bands.charcol(:) = color of fatbands
% bands.charscal = fatness of the fatbands = scaling factor 
% bands.plotbands = 0,1 (set to 0 if you don't want to see the normal bands
%                   but just the fatbands) 
%
% x = virtual x axis for band plots
% ax(:) = axis handles, one for each band

% extract data
npnts = size(bands.kpnt_pos,1);   % number of kpoints
nbands = size(bands.eval,2);    % number of bands

% generate a virtual x-axis:
% the unit of x correspond to real distances in the space of the k-point
% coordinates, which can be relative, Cartesian, etc. 
% as provided in the band data structure
%
% initialize variables: x(1) = 0 
x     = zeros(1,npnts);     
fat_x = zeros(1,2*npnts); 
fat_y = zeros(1,2*npnts); 

for kpnt = 2:npnts
    diffvec = bands.kpnt_pos(kpnt,:) - bands.kpnt_pos(kpnt-1,:);
    x(kpnt) = x(kpnt-1) + norm(diffvec);    
    
    % prepare x positions for plotting fatbands
    fat_x(kpnt) = x(kpnt);
    fat_x(2*npnts-kpnt+1) = x(kpnt);         
end

xmax = x(npnts);          % maximal x value
%aspratio = (bands.emax-bands.emin)/xmax;    % aspect ratio of the plot, needed to plot round circles

% now plot everything
hold on;

% if 'bands.bchar' is specified
% plot 'fatbands', indicating the band character,
% this is done by plotting one polygon for each band
if (~isscalar(bands.bchar))          
    for bnd = 1:nbands

        % original, unsmoothed data
        bchar_smooth = bands.bchar(:,bnd);
            
        % the plain band characters are quite noisy, therefore we smoothen
        % the band characters for each band along the complete k-path,
        % either by convolution or by spline; 
        % problem: for nice smoothing the data is strongly smeared out,
        % therefore one has to find the right balace between sufficient
        % smoothing and reduction of smearing; 
        % for optimized parameters both methods perform rather similar        
        
        % convolution = running averages
%          span = 8; % Size of the averaging window = number of data points
%          window = ones(span,1)/span;                          
%          bchar_smooth = convn(bands.bchar(:,bnd), window, 'same');                
        
        % cubic smoothing spline
        % smoothing parameter must be very close to
        % 1, otherwise the data is strongly smeared out
%         smoothing = 1-1E-5;
%         pp = csaps(x, bands.bchar(:,bnd), smoothing);
%         bchar_smooth = fnval(pp,x);
        
        % each polygon is defined by a serial list of edge points 
        % P(i) = [fat_x(i), fat_y(i)]
        % that are connected to define a filled 2D polygon; at each k-point
        % we generate the upper and the lower points of the polygon: 
        % fat_x = [x x_inverse]
        % fat_y = [eval-bchar (eval+bchar)_inverse]
        for kpnt = 1:npnts                               
            fat_y(kpnt) = bands.eval(kpnt,bnd) - bchar_smooth(kpnt)*bands.charscal;
            fat_y(2*npnts-kpnt+1) = bands.eval(kpnt,bnd) + bchar_smooth(kpnt)*bands.charscal;                        
        end
        
        % draw polygon
        fill(fat_x, fat_y, bands.charcol, 'FaceColor',bands.charcol, ...
            'EdgeColor','none', 'FaceAlpha', 0.5); %0.4
        
        
    end % bnds        
end

% plot black bands on top
if (bands.plotbands)
    if(isfield(bands, 'color')) 
        % this was implemented later, so make a check for backward
        % compatibility
        ax = plot(x, bands.eval, bands.linespec, 'Color', bands.color);
    else
        ax = plot(x, bands.eval, bands.linespec);
    end
else
    ax = 0;
end

% plot zeroline = EF
plot([0 xmax], [0 0], '--k');

% draw vertical separation lines
% the NUMBER of k-labels is used to determine the position of the
% separation lines because VASP always uses the SAME number of k-points
% along each k-line
nlabels = size(bands.klabels,2);
xind = [1 npnts/(nlabels-1)*(1:nlabels-1)];   % indices 'xind' of  x(xind) that correspond to the posn. of a k-point label
vtmp = ones(1,nlabels);                             % temporary row vector
xpnts = [x(xind); x(xind)];             % 2x(nlabels) dimensional matrix of x positions
ypnts = [bands.emin*vtmp; bands.emax*vtmp];         % 2x(nlabels) dimensional matrix of y positions
plot(xpnts, ypnts, '-k');  

% plot special point labels
elabel = bands.emin - (bands.emax-bands.emin)*.05;    % y position of the labels
text(x(xind), elabel*vtmp, bands.klabels, ...
    'HorizontalAlignment', 'center', 'Interpreter','latex');

% axes settings
axis([0 xmax bands.emin bands.emax]);   % visible range
set(gca, 'XTick', []);      % remove x axis ticks
ylabel('E-E_F (eV)');      % y label
box on                      % display boxed axes

hold off;
