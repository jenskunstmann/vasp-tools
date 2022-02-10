function [Eout, dos, idos] = el__CalculateDOS (Emat, weights, sigma)
% computes the gaussian smeared DOS from Eigenvalues Emat
% at points equally spaced by sigma in
% the range of Emat and padded to respect smearing of border values.
% Each Eigenvalue is smeared with a gaussian kernel with standard
% deviation sigma that is computed from Ei-10*sigma to E+10*sigma.
% weights contains a weight for each energy eigenvalue. For k-weighting,
% sum(weights(:)) must equal 1. For pDOS, multiply weights with weights, for
% example the absolute squares of eigenvectors.
%
% Emat(n,k)    = eigenvalues in some energy unit, n=band index, k=k-point index 
% weights(n,k) = corresponding weights
%                weights(1,:) for normal dos = only k-point weights
%                weights(:,:) for projected dos, where the n-related weight
%                correspond to the orbital character
% sigma        = Gaussian smearing width in the same energy unit 
%                
%
% if correct k-point weights are provided:
% [dos(:)]  = states/(energy unit)/cell
% [Eout(:)] = (energy unit)
% [idos(:)] = states/cell
%
% This is the result of the function:
% k=exp(-((E-E_i)/sigma)^2)
% dos = sum_i C_i * k_i / sum(k_i) / sigma 
%
% For projected DOS, pass weights as a matrix of size(Emat) and appropriate
% weights, e.g.
%   weights = |<band|orbital>|^2 = |EV(k,band)|^2
%   >> weights = norm_EV(EV);
%   >> weights = permute(weights(orbital,:,:,spin),[2 3 1 4]).^2;
%
% If k-Points do not have equal weights, multiply columns of weights with those
% weights (set weights=1 if you don't want to specify it):
%   >> CW = weights .* W(ones(1,size(weights,1)),:)
%
% Copyright 2010, Felix Zörgiebel


if nargin==0, disp('Usage:[Eout, dos, idos] = dos_smear (Emat, weights, sigma, Eaxis)'); return, end

% sigmafault weights: 1/N
if nargin<2 || isempty(weights)
    disp('Assuming constant weights for all k-points!')
    weights=1/size(Emat,2); % #k-values=#unit-cells=size(Emat,2)
end

if size(weights,1)==1
    weights=repmat(weights,size(Emat));
end

mi=min(Emat(:));
ma=max(Emat(:));

if nargin<3 || isempty(sigma)
    sigma = max((ma-mi)/(round(numel(Emat)/1.3)-1),(ma-mi)/799)
end

%if nargin<4 || isempty(ne), ne=size(Emat,1); end

w=10; % this is the width of the actually computed part of the kernel function. It's more than sufficient (gauss goes below 1e-6)

Eout = (floor(mi/sigma)-w:ceil(ma/sigma)+w)*sigma;

% kernel function
kernel=@(x) exp(-(x/(sigma)).^2);

% sum up kernels for all eigenvalues
dos = zeros(size(Eout));
for i=1:size(Emat,1)
    for j=1:size(Emat,2)
        E=Emat(i,j);
        idx=(0:2*w)+round((E-mi)/sigma)+1;
        k=kernel(Eout(idx)-E);
        k=k/sum(k); % norm kernel to one state
        dos(idx)=dos(idx)+k*weights(i,j);
    end
end

% scale dos to states/(energy unit)/cell; 
% factor 2 = spin
dos=2*dos/sigma; 

idos = cumsum(dos)*sigma;

end