function [omega, eps1, eps2, n, kappa, alpha, loss] = vasp__getOpticalCoefficients(varargin)
% read in the vasprun.xml or dielectric.mat, extract the dielectric function, calculate the
% optical coefficients and return them 
%
% omega(npnts,1)   .... frequencies in eV
% eps1(npnts, comp) ... real part of the dielectric function
% eps2(npnts, comp) ... imaginary part of the dielectric function
% n(npnts, comp)    ... refractive index
% kappa(npnts, comp)... extrinction coefficient
% alpha(npnts, comp)... absorption coeff
% loss(npnts, comp) ... electron energy loss function
% 
%           comp = 1 = xx
%           comp = 2 = yy
%           comp = 3 = zz
%           comp = 4 = xy
%           comp = 5 = yz
%           comp = 6 = zx   

global SYS


%%%%%%%%%% OPTIONS
epsdir = 'dos';     % subdirectory that contains the DF
%%%%%%%%%%%

const = PhysicalConstants();

% read vasprun.xml file
% eps2(npnts, comp) : eps2inary part of the dielectric function
% eps1(npnts, comp) : eps1 part of the dielectric function, with 'npnts'
%                        energy points, and comp = 1..7
%                           1 = energy
%                           2 = xx
%                           3 = yy
%                           4 = zz
%                           5 = xy
%                           6 = yz
%                           7 = zx    

% define file source by arguments
switch(nargin)
    case(1)
        matfile = varargin{1}      
        
    otherwise
        % reading in the data that was extracted before with the python tool:
        % ~/prog/python/vasp-scripts/vgetdielectric.py  
        matfile = sprintf('%s/%s/%s', SYS.path, epsdir, 'dielectric.mat')             
end

% read in data
load(matfile);   

% reading in the xml directly (very slow for big files)
%xmlfile = sprintf('%s/%s/%s', SYS.path, epsdir, 'vasprun.xml')
%[real, imag] = vasp__readDielectricFunction(xmlfile);



% extract energy, and the relevant component of the DF
omega = real(:,1);
eps1  = real(:,2:7);
eps2  = imag(:,2:7);

n = sqrt(eps1 + sqrt(eps1.^2 + eps2.^2))/sqrt(2);
kappa = sqrt(-eps1 + sqrt(eps1.^2 + eps2.^2))/sqrt(2);


%R  = ((n-1).^2 + kappa.^2) ./ ((n+1).^2 + kappa.^2);

loss = eps2 ./ (eps1.^2 + eps2.^2);

c = 4*pi*const.e/const.h/const.c/100;
alpha = c*kappa.*omega(:,ones(1,6));  % [alpha] = 1/cm




