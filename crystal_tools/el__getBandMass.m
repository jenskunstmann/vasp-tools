function [m] = el__getBandMass(x, bands, band_id, use_npnts)
% calculate the effective/band masses at the special points along one band
%
% x(:)      = virtual x axis for band plots; MAKE SURE THIS HAS THE UNIT [1/Ang] !!!
%
% bands.emax = plotting range
% bands.emin 
% bands.klabels{} =  special point labels
% bands.eval(kpnt,band) = band energies with Fermi level at E=0
% bands.kpnt_pos(kpnt,:) = positions of kpoints (in any given unit)
% bands.linespec = 'linespec' of the lines of the bands
% bands.bchar(kpnt,band) = band characters = fatbands
% bands.charcol(:) = color of fatbands
% bands.charscal = fatness of the fatbands = scaling factor 
% bands.plotbands = 0,1 (set to 0 if you don't want to see the normal bands
%                   but just the fatbands) 
%
% band_id   = band number used for the determination of the band masses
% use_npnts = number of points to use for the harmonic fit
%
% m(:) = vector with effective masses in units of [m_e]
%
% THERE IS STILL AN OPEN ISSUE WITH THE RIGHT UNIT: 
% its unclear why to divide the masses by 2
% However for MoS2 mono, bilayers, GaAs and Si the results were in very
% good agreement with literature

% physical constants
c = PhysicalConstants();

% extract data
npnts = size(bands.kpnt_pos,1);     % number of kpoints
nbands = size(bands.eval,2);        % number of bands
nlines = size(bands.klabels,2) - 1; % number of k-lines
pntsperline = npnts/nlines;         % points per k-line

hold on

% determine band masses at the two boudary points
for line = 1:nlines
    
    % select one k-line
    indx_kline_start = (line-1)*pntsperline + 1;
    indx_kline_end   = (line)*pntsperline;    
    %x(indx_kline_start:indx_kline_end)
    
    % pick left and right parts of the k-line
    indx_left_start  = indx_kline_start;
    indx_left_end    = indx_kline_start + use_npnts - 1;
    indx_right_start = indx_kline_end - use_npnts + 1;
    indx_right_end   = indx_kline_end;
    
    % get the distances and energies
    k_left  = x(indx_left_start:indx_left_end);
    k_right = x(indx_right_start:indx_right_end);
    eval_left  = bands.eval(indx_left_start:indx_left_end,   band_id)';    
    eval_right = bands.eval(indx_right_start:indx_right_end, band_id)';    
    
    % display values for other use
    xy_left = [k_left' eval_left']
    xy_right = [k_right' eval_right']

    % LEFT point -> k-dispersion towards right
    % fit parabola: y = c*(x-x0)^2 + b
    % see generalfit() for explanation of the technique    
    x_m_x0_sq = (k_left - k_left(1)).^2;                   % = (x-x0)^2
    y_m_b = eval_left - eval_left(1)*ones(1, use_npnts);   % = y-b    
    const = y_m_b / x_m_x0_sq;                             % least square fit     
    
    % make harmonic fit: coefficients in order of decending power        
    %coeff = polyfit(k_left, eval_left, 2);
    
    % WHY DO WE NEED TO DIVIDE BY 2? -> (d2y/dk2) = 2*c 
    A = (c.h/2/pi)^2*1E20/(c.me * c.e)/2;  % = 7.62/2  
    m(1,line) = A/const;
    %m(1,line) = A/coeff(1);

    % plot the current fit
    y = const * x_m_x0_sq + eval_left(1);    % evaluate fitted function
    %y = polyval(coeff, k_left);
    plot(k_left, y, '-r' )

    
    % RIGHT point -> k-dispersion towards left
    % fit parabola: y = ax^2 + b
    % see generalfit() for explanation of the technique   
    x_m_x0_sq = (k_right - k_right(use_npnts)).^2;                  % = (x-x0)^2
    y_m_b = eval_right - eval_right(use_npnts)*ones(1, use_npnts);   % = y-b    
    const = y_m_b / x_m_x0_sq;                                      % least square fit
    y = const * x_m_x0_sq + eval_right(use_npnts);    % evaluate fitted function
    plot(k_right, y, '-b')    
    m(2,line) =  A/const; 
              
%     % make harmonic fit: coefficients in order of decending power
%     coeff = polyfit(k_right, eval_right, 2);
%     m(2,line) =  A/coeff(1);        
%     
%     % plot the current fit
%     %y = polyval(coeff, k_right);
%     plot(k_right, y, '-b' )              
    
end

display(sprintf('\nband masses along band nr=%d, ordered as marked in the plot, unit [m_e]', band_id))
% output
m = reshape(m, nlines*2, 1)';


