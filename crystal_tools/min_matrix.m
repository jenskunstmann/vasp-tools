function [minimum, line, column] = min_matrix(A)
% determine the minimum value of a matrix and the corresponding matrix
% indices

[minvec ,lineindx]    = min(A);
[minimum, column] = min(minvec);
line = lineindx(column);
