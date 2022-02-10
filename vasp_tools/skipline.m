function skipline(fh, n)
% skip n lines when reading a file specified by the file handle 'fh'
% empty space, tab and newline are totally ignored

for i = 1:n
    head = fgetl(fh);
    %disp(head)
end;