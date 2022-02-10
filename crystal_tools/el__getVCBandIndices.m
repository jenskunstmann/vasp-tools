function [vbandindx, cbandindx, mbandindx] = el__getVCBandIndices(bands)
% determine the indices of valence, conduction and metallic bands and the
% maxima and minia of each band
% 
% bands = structure for plotting bands, with EF = 0 eV 
%
% vbandindx(:)
% cbandindx(:)
% mbandindx(:) = array with the indices of all valence/conduction/metallic
%                bands

% prepare variables
kpnt_pos = bands.kpnt_pos;
eval = bands.eval;
[nkpnts, nbands] = size(eval);


%bndmax_x = zeros(1,nbands); bndmax_y = zeros(1,nbands);
%bndmin_x = zeros(1,nbands); bndmin_y = zeros(1,nbands);
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

% get the VBM, CBM and the band gap
vbandmax_y = bndmax_y(vbandindx(vbandcnt));  % VBM
%vbandmax_x = bndmax_x(vbandindx(vbandcnt));  
cbandmin_y = bndmin_y(cbandindx(1));         % CBM 
%cbandmin_x = bndmin_x(cbandindx(1));          
bandgap = cbandmin_y - vbandmax_y

