function [minimum, line, column] = max_matrix(A)
% determine the minimum value of a matrix and the corresponding matrix
% indices

[minvec ,lineindx] = max(A);
[minimum, column]  = max(minvec);
line = lineindx(column);
