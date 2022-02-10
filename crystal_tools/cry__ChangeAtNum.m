function newpatch = cry__ChangeAtNum(oldpatch, AtNumShift)
% shift the atomic numbers of all atoms by 'atnumshift';
% the purpose of this is to change their appearance (color/diameter) when
% plotting them, because PlotAtomicStructure() assigns a different
% color/diameter to every element
% MAKE SURE: that the atomic number does not exced 120

newpatch = oldpatch;
newpatch.atomnum = oldpatch.atomnum + AtNumShift;
end