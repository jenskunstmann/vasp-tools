function [nwannier] = vasp__readWannier90Wout(filename)
% read wannier90.wout file and return some values 
%
% syntax: [nwannier] = vasp__readWannier90Wout(filename)
% 
% nwannier = number of Wannier orbitals


% read file content into string
contstr = fileread(filename);
contstrlen = length(contstr);

% get some general parameters
nwannier = getvalue(contstr, 'Number of Wannier Functions               :');


function value = getvalue(str, searchstr)
% return numerical value after keyword

index = strfind(str, searchstr);
if isempty(index)
    error(sprintf('keyword <%s> not found.', searchstr))
end

startindex = index(1) + length(searchstr);
endindex   = startindex + 50;
[value, count] = sscanf(str(startindex:endindex), ' %f ');
if count ~= 1
    error('error reading file');         
end;