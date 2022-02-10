function [estart eend] = el__getBandGaps(e, dos)
% extract all band gaps from a density of states (DOS) 
%
% e(1..npnt) : column vector with the energies, E_F = 0
% dos(1..npnt) : 1D column vector of the dos
%
% estart(:) : energies where a band gap starts
% eend(:)   : energies where a band gap ends
% 
% the actual gaps are simply: eend - estart

lendos = length(dos);
emin = min(e);
emax = max(e);

% find potential gaps
isgap = false;
gapcount = 0;
for index = 1:lendos
    % beginning of gap
    if (isgap == false) && (dos(index) == 0)
        isgap = true;
        gapcount = gapcount + 1;
        estart(gapcount) = e(index);
        eend(gapcount) = NaN;
        gap(gapcount) = NaN;
        
    % end of gap    
    elseif (isgap == true) && (dos(index) ~= 0)   
        eend(gapcount) = e(index-1);
        gap(gapcount) = e(index-1) - estart(gapcount);        
        isgap = false;               
    end
end

% remove fake gaps at the boundaries, zero gaps, and NaN elements
if estart(1) == emin
    gap(1)    = NaN;
end 

% mask out the unwanted elements, i.e. NaN or zero gaps
delcount = 0;
for gapindx = 1:gapcount
    if (gap(gapindx) == 0) || (isnan(gap(gapindx)))
        delcount = delcount + 1;
        delindx(delcount) = gapindx;
    end
end

% remove unwanted elements
estart(delindx) = [];
eend(delindx)   = [];
gap(delindx)    = [];
        
% potential output
%[estart'  eend' gap']

