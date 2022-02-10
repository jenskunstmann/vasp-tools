function smoothdata = RunningAverage(plaindata, swidth, periodic)
% smoothens a curve via a running average
%
% plaindata(:) = vector with raw data, 
%                must be a line vector NOT column vector
% swidth = width of average = number of data points
% periodic = 1: the raw data is periodic in the given range (the
%               periodicity must match), the running average is complete
%               throughout the whole data range (due to padding)
%          = 0: the raw data is not periodic and the running average is
%               incomplete near the boundaries
%
% smoothdata(:) = smoothed data

npnts = length(plaindata);

% convolution with a constant and normalized window function = running average
% define convolution window
window = ones(1,swidth)/swidth;

%wnorm = sum(sum(window));
%if wnorm ~= 1, display('window function must be normalized to one to be norm conserving'), end

if periodic
    % do padding (extend the scale by the width of the window) to get the
    % convolution in the central part right 
    plaindata_pad = [plaindata(npnts-swidth+1:npnts) plaindata plaindata(1:swidth)];

    % do the actual convolution
    smoothdata_pad = convn(plaindata_pad, window, 'same');

    % unpad the convoluted curve
    smoothdata = smoothdata_pad(swidth+1:npnts+swidth);    
else
    smoothdata = convn(plaindata, window, 'same');
end


