function el__plotBandNumbers(bands, virtax, color, lshift)
% this routine can be added within a band plot to add band numbers
% 
% bands ... bands structure
% virtax ... virtual x axis
% color(:) ... rgb color vector
% lshift(x,y) ... shift of text relative to the original position

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
            % plot the labels outside of the actual horizontal plot
            if kpnt == 1
                text(double(virtax(kindex)-lshift(1)), double(energy+lshift(2)), ...
                     int2str(n), 'Color', color, 'HorizontalAlignment', 'right');              
            else
                text(double(virtax(kindex)+lshift(1)), double(energy+lshift(2)), ...
                     int2str(n), 'Color', color, 'HorizontalAlignment', 'left');           
            end            
        end
    end
end