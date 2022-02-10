function transitions = el__getBandTransitions(bands, x)
% determine all direct and indirect transitions between the valence and the
% conduction band, plot the transitions points on top of the band structure
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
% kpnt_pos(kpoint, component): position of the k-point in reduced
%                              coordinates,  i.e. in units of the 
%                              reciprocal lattice vectors
%
% ROUTINE IS PROBABLY NOT WORKING BECAUSE el__getVCBandIndices() was cut
% out.



%%%%% USER DATA %%%%%%%%
text_shift = [0.01 -.3];
%%%%%%%%%%%%%%%%%%%%%%%%%5

% prepare variables
kpnt_pos = bands.kpnt_pos;
eval = double(bands.eval);  % some routines cannot cope with single precision numbers
[nkpnts, nbands] = size(eval);


%bndmax_x = zeros(1,nbands); bndmax_y = zeros(1,nbands);
%bndmin_x = zeros(1,nbands); bndmin_y = zeros(1,nbands);
vbandindx = []; vbandcnt = 0;
cbandindx = []; cbandcnt = 0;
mbandindx = []; mbandcnt = 0;
for ibnd = 1:nbands
    [bndmax_y(ibnd), indx] = max(eval(:,ibnd));
    bndmax_x(ibnd) = x(indx);
    [bndmin_y(ibnd), indx] = min(eval(:,ibnd));   
    bndmin_x(ibnd) = x(indx);

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

% get the VBM, CBM and the band gap
vbandmax_y = bndmax_y(vbandindx(vbandcnt));  % VBM
vbandmax_x = bndmax_x(vbandindx(vbandcnt));  
cbandmin_y = bndmin_y(cbandindx(1));         % CBM 
cbandmin_x = bndmin_x(cbandindx(1));          
bandgap = cbandmin_y - vbandmax_y

% find the min/max values for direct/indirect transitions at the band gap
[v_xm, v_ym, v_kindx] = gen__FindPeaks(x, eval(:, vbandindx(vbandcnt)), 'max');
[c_xm, c_ym, c_kindx] = gen__FindPeaks(x, eval(:, cbandindx(1)), 'min');

% add the VBM and CBM in case they are not yet included (e.g. if they are
% at the boundaries of the plot) and remove double occurances (unique(), the option
% 'stable' suppresses sorting)
v_xm = unique([v_xm;  vbandmax_x], 'stable'); v_ym = unique([v_ym;  vbandmax_y], 'stable');
c_xm = unique([c_xm;  cbandmin_x], 'stable'); c_ym = unique([c_ym;  cbandmin_y], 'stable');
%v_xm = unique([v_xm;  vbandmax_x]); v_ym = unique([v_ym;  vbandmax_y]);
%c_xm = unique([c_xm;  cbandmin_x]); c_ym = unique([c_ym;  cbandmin_y]);

% generate text labels
for i = 1:length(v_xm) 
    v_text{i} = num2str(i);
end
for i = 1:length(c_xm) 
    c_text{i} = num2str(i);
end

% plot transitions points and text labels
hold on
plot(v_xm, v_ym, 'o', 'MarkerSize', 10);
text(v_xm + text_shift(1), v_ym + text_shift(2), v_text, 'Color', [0 0 1])
plot(c_xm, c_ym, 'ro', 'MarkerSize', 10);
text(c_xm + text_shift(1), c_ym + text_shift(2), c_text, 'Color', [1 0 0])

% print k-position of critical (min/max energy) k-points in units of the
% reciprocal unit cell vectors (boundary points are not included)
display(sprintf('k-positions of transition points in units of the reciprocal unit cell vectors:'))
kpos_valence    = kpnt_pos(v_kindx,:)
kpos_conduction = kpnt_pos(c_kindx,:)

% calculate the transition energies
transitions = zeros(length(c_ym), length(v_ym));
for vi = 1:length(v_ym)
    for ci = 1:length(c_ym)
        transitions(ci,vi) = c_ym(ci) - v_ym(vi);
    end
end

display(sprintf('1 2 3 .... -> points in valence band\n2\n3\n:\n:\n|\nv points in conduction band'))
transitions

