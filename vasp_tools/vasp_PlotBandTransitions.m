function vasp_PlotBandTransitions()
% make a plain bandstructure plot for VASP

global SYS

%%%%%% USER PART %%%%%%%%%%
% energy window for the transition energies and the spectral function
emin = 3; 
emax = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read EIGENVAL file : 
file_eigenval = sprintf('%s/%s',SYS.path,SYS.eigenval)
[kpnt_pos eval_up eval_down] = vasp__readEIGENVAL(file_eigenval);
% read PROCAR file
% file_procar = sprintf('%s/%s',SYS.path,SYS.procar)
% [kpnt_pos, eval_up, bandchar] = vasp__readPROCAR(file_procar, SYS.lorbit);
% eval_down = [1];

% read imaginaly part of dielectric function
%subpath = 'dos/imag.dat';
%file_imag = sprintf('%s/%s', SYS.path, subpath)
%I = importdata(file_imag);
[omega, eps1, eps2, n, kappa, alpha, loss] = vasp__getOpticalCoefficients();

omega = omega(:,1); 
eps2 = eps2(:,1);
% find the peak positions
%[pks,locs] = findpeaks(eps2);
%peaks = [omega(locs) pks];


% get Fermi energy from DOSCAR
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar);
efermi = vasp__getEFermi(file_doscar);

% compile the relevant data for plotting
bands.klabels = SYS.klabels;          % special point labels
bands.kpnt_pos = kpnt_pos;        % positions of the kpoints  
bands.bchar = 0;                  % no band characters
bands.charscal = 0;               % no fatness for the band characters = fatbands
bands.charcol = 0;                % color of the fatbands
bands.plotbands = 1;

% plot bands 
if(size(eval_down) == [1 1])     % non-spin-polarized if 'eval_down' is empty
    bands.eval = eval_up - efermi;     
    bands.linespec = '-k';   
    
    % plot normal band structure with band numbers
    subplot(1,3,1)
    bands.emin = -2; %SYS.emin;
    bands.emax = 4; %SYS.emax;                % plotting range
    x = vasp__plotBands(bands);  
    vasp__plotBandNumbers(bands, x, [0 0 1], [0.02 0]); 
    title(SYS.ID,'Interpreter','none');
    
    % plot band transitions 
    subplot(1,3,2)
    bands.emin = emin;    
    bands.emax = emax;       % plotting range
    [tbands tlabels] = getBandTransitions(bands);
    x = vasp__plotBands(tbands); 
    PlotBandTransitions(tlabels, tbands, x, [0 0 1]);
    hold on
    % plot peak positions as horizontal lines
    xpos(1,:) = zeros(length(pks),1);
    xpos(2,:) = max(x)*ones(length(pks),1);
    plot(xpos, [peaks(:,1) peaks(:,1)]' , '-r')
    title('Band Transitions','Interpreter','none');
    ylabel('E_c - E_v (eV)')
    
    % plot absorption
    subplot(1,3,3)
    hold on
    plot(pks, omega(locs), 'or');
    plot(eps2, omega, 'k');  % plain
    title('Dielectric Function');
    xlabel('\epsilon_2')
    ylabel('Energy (eV)')  
    ylim([emin emax])
    box on
    
else                             % spin-polarized
    % spin-up bands
    bands.eval = eval_up - efermi;     
    bands.linespec = '-r';             
    tbands = getBandTransitions(bands);
    vasp__plotBands(tbands);

    % spin-down bands
    bands.eval = eval_down - efermi;   
    bands.linespec = '-b';             
    tbands = getBandTransitions(bands);
    vasp__plotBands(tbands);
end;

% save band energies in an ASCII file
%ematrix = eval_up - efermi;
%save('bands.dat', 'ematrix', '-ASCII')


function [tbands tlabels] = getBandTransitions(bands)
% determine all vertical transitions between the valence and the
% conduction bands and return as 'bands' structure
% 
% bands.emax = plotting range
% bands.emin 
% bands.klabels{} =  special point labels
% bands.eval(kpnt,band) = band energies with Fermi level at E=0
% bands.kpnt_pos(kpnt,:) = positions of kpoints
% bands.linespec = 'linespec' of the lines of the bands
% bands.bchar(kpnt,band) = band characters = fatbands
% bands.charcol(:) = color of fatbands
% bands.charscal = fatness of the fatbands = scaling factor 
% bands.plotbands = 0,1 (set to 0 if you don't want to see the normal bands
%                   but just the fatbands) 
% x = virtual x axis of the bands plot

% initialize the output
tbands = bands;

eval = bands.eval;
[nkpnts, nbands] = size(eval);

% determine the indices of valence, conduction and metallic bands and the
% maxima and minia of each band
%bndmax_x = zeros(1,nbands); 
bndmax_y = zeros(1,nbands);
%bndmin_x = zeros(1,nbands); 
bndmin_y = zeros(1,nbands);
vbandindx = []; vbandcnt = 0;
cbandindx = []; cbandcnt = 0;
mbandindx = []; mbandcnt = 0;
for ibnd = 1:nbands
    [bndmax_y(ibnd), indx] = max(eval(:,ibnd));
    %bndmax_x(ibnd) = x(indx);
    [bndmin_y(ibnd), indx] = min(eval(:,ibnd));   
    %bndmin_x(ibnd) = x(indx);

    % valence bands
    if (bndmax_y(ibnd) < 0) && (bndmin_y(ibnd) < 0)
        vbandcnt = vbandcnt + 1;
        vbandindx(vbandcnt) = ibnd;
        
    % conduction bands 
    elseif (bndmax_y(ibnd) > 0) && (bndmin_y(ibnd) > 0)
        cbandcnt = cbandcnt + 1;
        cbandindx(cbandcnt) = ibnd;
        
    % bands cutting EF    
    else
        mbandcnt = mbandcnt + 1;
        mbandindx(mbandcnt) = ibnd;
        
    end

end

% initialze
tval   = zeros(nkpnts, vbandcnt*cbandcnt);      % the actual transition energies
tlabels = zeros(nkpnts, vbandcnt*cbandcnt, 2);   % a 'label' that tells us what bands where considered
bandcnt = 0;
% calcualte all transition energies
for vbi = 1:vbandcnt
    for cbi = 1:cbandcnt    
        bandcnt = bandcnt + 1;
        for kpnt = 1:nkpnts            
            tval(kpnt,bandcnt) =  eval(kpnt,cbandindx(cbi)) - eval(kpnt,vbandindx(vbi));
            tlabels(kpnt,bandcnt,:) = [cbandindx(cbi), vbandindx(vbi)];
        end
    end
end

tbands.eval = tval;


function PlotBandTransitions(tlabels, bands, virtax, color)
% this routine plots the band indices of the transitions
% 
% bands ... bands structure
% virtax ... virtual x axis
% color(:) ... rgb color vector

% label shift
lshift = [0.02 0];

% extract data
nkpnts = size(bands.kpnt_pos,1);    % number of kpoints
nbands = size(bands.eval,2);        % number of bands
nlabels = size(bands.klabels,2);    % number of special points

% plot band numbers
xind = [1 nkpnts/(nlabels-1)*(1:nlabels-1)];  % indices 'xind' of  x(xind) that correspond to the posn. of a k-point label
for kpnt = 1:nlabels
    for n = 1:nbands
        kindex = xind(kpnt);
        energy = bands.eval(kindex,n);

        % do not plot labels outsinde of the vertical plotting range = ugly
        if (energy > bands.emin) && (energy < bands.emax)   
            
            % plot the labels outside of the actual plot
            if kpnt == 1
                text(virtax(kindex)-lshift(1), bands.eval(kindex,n), ...
                     sprintf('%d-%d',tlabels(kindex,n,:)), 'Color', color, ...
                     'HorizontalAlignment', 'right');              
            else
                text(virtax(kindex)+lshift(1), bands.eval(kindex,n), ...
                     sprintf('%d-%d',tlabels(kindex,n,:)), 'Color', color, ...
                    'HorizontalAlignment', 'left');                
            end
        end
    end
end


