function [xm, ym, xindx] = gen__FindPeaks(x, y, type)
% depending on 'type' find the up-peaks or down-peaks of the curve given by x,y
%
% xm, ym ... the x,y coordinates of the minima/maxima
% xindx ... index along x == k point index
%
% uses findpeaks() from signal processing tool box.

switch(type)
    case 'min'
        [ym,xindx] = findpeaks(1-y);
        ym = 1-ym;
        xm = x(xindx)';
        
    case 'max'
        [ym,xindx] = findpeaks(y);
        xm = x(xindx)';        
        
    otherwise
        error('Type of operation not defined. Use type = <min> or <max>.')
end
