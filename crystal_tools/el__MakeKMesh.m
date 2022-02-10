 function [kpoints, weights] = el__MakeKMesh(subdiv, shift)
% generate the kpoints of a equally spaced k-mesh with origin at Gamma, if
% shift = [0 0 0]  
%         
% subdiv = [n1 n2 n3] 
%        = number of subdivisions along each rlv
% shift  = [d1 d2 d3] 
%        = shift of the k-meh off the center of the BZ in units of the rlv
%          meaningful range: 0 < di < 1
%
% kpoints(1..nkpnts, comp) = coordinates of the k-points in units of the rlv
% weights(1..nkpnts)       = k-point weights, 
%                            currently each k-point has the same weight
%                            sum(weights) == 1 (normalization condition)

N1 = subdiv(1);
N2 = subdiv(2);
N3 = subdiv(3);
nkpnts = N1*N2*N3;  % total number of k-points

% equally spaced mesh with origin at Gamma, if shift = [0 0 0] 
% = 'Gamma' option in VASP KPOINTS
cnt = 0;
kpoints = zeros(nkpnts, 3);
for i1 = 0:N1-1
    for i2 = 0:N2-1
        for i3 = 0:N3-1
            cnt = cnt +1;
            % position in units of the rlv, mapped back to the 1st cell
            kpoints(cnt, :) = mod([i1/N1 + shift(1), i2/N2 + shift(2), i3/N3 + shift(3)],1);
        end 
    end 
end

% make k-point weights, each point is weighted as 1/nkpnts
weights = ones(nkpnts,1)/nkpnts;
