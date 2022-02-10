function vasp_PlotOpicalCoeffs_all()
% plots the frequency dependent dielectric function (DF) and other related
% optical coefficients

global SYS


%%%%%%%%%% OPTIONS
comp = 1;           %
%           comp = 1 = xx
%           comp = 2 = yy
%           comp = 3 = zz
%           comp = 4 = xy
%           comp = 5 = yz
%           comp = 6 = zx   
window = 0.3;       % smoothing window in eV, must be > 0
%%%%%%%%%%%

%matfile = '/exports-hal/users/home/jkunstmann/res/MoS2/DoubleLayers/13.2deg/initstruc-Mo-S-new/relax-algo2/dos-Kpnt/dielectric.mat'
%matfile = '/home/jk/res/TMD/HS/0.05%/21deg/relax-algo2/dos-Gpnt/dielectric.mat'
[omega, eps1, eps2, n, kappa, alpha, loss] = vasp__getOpticalCoefficients();

% extract energy, and the relevant component of the DF
omega = omega(:,comp);
eps1  = eps1(:,comp);
eps2  = eps2(:,comp);
n     = n(:,comp);
kappa = kappa(:,comp);
alpha = alpha(:,comp);
loss  = loss(:,comp);

% extract data
npnts = length(eps2);
domega = omega(2)-omega(1);

% sunbplot scheme
nlines = 4;
ncolumns = 2;


%%%%%%% eps1  part of DF
subplot(nlines, ncolumns ,1)
hold on
plot(omega, eps1, 'r');  
ylabel('\epsilon_1')
title('dielectric function, real part ')
plotoptions();

%%%%%%% imaginary part of DF
subplot(nlines, ncolumns ,2)

% window_pnts = ceil(window/domega) % number of data points the window corresponds to
% eps2_av = RunningAverage(eps2, window_pnts, 0);
% 
% % check norms
% sum(eps2/npnts) % norm of the bare data
% sum(eps2_av/npnts) % norm of the smoothed data

% find the peaks
%[pks,locs] = findpeaks(eps2);
%peaks = [omega(locs) pks];
%peaks(:,1)
% [pks,locs] = findpeaks(eps2_av);
% peaksav = [omega(locs); pks]'
% 
% plot absorption
hold on
plot(omega, eps2, 'r');  % plain
%plot(peaks(:,1), peaks(:,2), 'or');
% plot(omega, eps2_av, '--k');  % snmoothed
% plot(peaksav(:,1), peaksav(:,2), 'ok');
ylabel('\epsilon_2')
title('dielectric function, imag. part')
plotoptions();


%%%%%%% diffractive index
subplot(nlines, ncolumns ,3)
%n = sqrt(eps1 + sqrt(eps1.^2 + eps2.^2))/sqrt(2);

hold on
plot(omega, n, 'r');  
ylabel('n')
xlabel('Energy (eV)')
title('refractive index')
plotoptions();

%%%%%%% extinction coefficient
subplot(nlines, ncolumns ,4)
%kappa = sqrt(-eps1 + sqrt(eps1.^2 + eps2.^2))/sqrt(2);

hold on
plot(omega, kappa, 'r');  
ylabel('\kappa')
title('extinction coefficient')
plotoptions();

%%%%%%% absorption coefficient
const = PhysicalConstants();
subplot(nlines, ncolumns ,6)
%c = 4*pi*const.e/const.h/const.c/100;
%alpha = c*kappa.*omega;  % [alpha] = 1/cm

hold on
plot(omega, alpha, 'r');  
ylabel('\alpha (cm^{-1})')
title('absorption coefficient')
plotoptions();

%[pks,locs] = findpeaks(alpha);
%peaks_alpha = [omega(locs) pks];
%peaks_alpha(:,1)
%plot(peaks_alpha(:,1), peaks_alpha(:,2), 'or');


% monolayer fractional change in reflectance
% no real difference in the C peal frequency found when considering dRbyR
% or only alpha
% see Mak, K. F., Sfeir, M. Y., Wu, Y., Lui, C. H., Misewich, J. a., &
% Heinz, T. F. (2008). Measurement of the Optical Conductivity of Graphene.
% Physical Review Letters, 101(19), 196405. doi:10.1103/PhysRevLett.101.196405 
subplot(4,2,8)
% linear approximation for the frequency dep. of the refective index of
% glass 
% between 2.4 and 3.3 eV (C peak), see Mark Fox, opti. prop of solids, p. 46
x=[2.4814 3.3085]; y = [1.462 1.472]
ns = (y(2)-y(1))/(x(2)-x(1))*(omega - x(1)) + y(1);
dRbyR = alpha*4./(ns.^2 - 1);

% peaks
hold on
%[pks,locs] = findpeaks(dRbyR);
%peaks_dRbyR = [omega(locs) pks];
%peaks_dRbyR(:,1)
%plot(peaks_dRbyR(:,1), peaks_dRbyR(:,2), 'or');

% plot
plot(omega, dRbyR, 'r');  
ylabel('\Delta R/R (a.u.)')
plotoptions();

% %%%%%%% reflectivity, according to Fresnel equations for normal incidence
% subplot(4,2,8)
% R  = ((n-1).^2 + kappa.^2) ./ ((n+1).^2 + kappa.^2);
% 
% hold on
% plot(omega, R, 'r');  
% ylabel('R')
% title('reflectivity')
% xlabel('Energy (eV)')
% plotoptions();
% 
% %%%%%%% EELS spectrum
% subplot(4,2,7)
% % eps2/(eps1^2 + eps2^2)
% %loss = eps2 ./ (eps1.^2 + eps2.^2);
% 
% hold on
% plot(omega, loss, 'r');  
% title('EELS energy loss function')
xlabel('Energy (eV)')
% plotoptions();

% big title
axes('Position',[0 0 1 1],'Visible','off');
text(0.45, 0.97, SYS.ID)

% peaks 
%[peaks(1:2) peaks_alpha(1:2) peaks_dRbyR(1:2)]


function plotoptions()
% options
xlim([1 4])
box on

