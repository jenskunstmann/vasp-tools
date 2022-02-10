function newpatch = cry__ShiftAtom(oldpatch, atomIDs, shift)
% shift the atom(s) with ID 'atomIDs(:)' by the vector 'shift'

newpatch = oldpatch;
nshatoms = length(atomIDs);
newpatch.atompos(atomIDs,:) = newpatch.atompos(atomIDs,:) + shift(ones(nshatoms,1),:);

end

% for at = 1:length(atomIDs)
%     newpatch.atompos(atomIDs(at),:) = newpatch.atompos(atomIDs(at),:) + shift;
% end