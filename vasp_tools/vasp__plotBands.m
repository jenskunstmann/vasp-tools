function [x, ax] = vasp__plotBands(bands)
% simple wrapper, so we do not have to change all the routines

[x, ax] = el__plotBands(bands);
