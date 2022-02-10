function[iangle]=b_myangle(av1,av2)
%b_myangle calculates the angle between two vectors av1 and av2
%USAGE: [iangle]=b_myangle(av1,av2)
%
if (size(av1)==size(av2))
    iangle = acosd(dot(av1,av2) / (norm(av1) * norm(av2)));
else
    error('vectors have to have the same dimension')
end
%END