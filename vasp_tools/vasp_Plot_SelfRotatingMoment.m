function [L] = vasp_Plot_SelfRotatingMoment()
% Calculates self-rotating angular momentum (Lself) from matrix elements 
% and band energies read from VASP's WAVEDERF 
% Assumes that WAVEDERF contains only one k point
% plots the convergence of L with respect to the number of bands included
% in the summation

global SYS

% read in files and data
file = 'WAVEDERF_K+';   % potentially add this file location to the SYS variable
filename = sprintf('%s/%s', SYS.path, file)
[matrix_elements,band_energies,occupations] = vasp__readWAVEDERF(filename);

% extract dimensions and VBM
nbands   = length(band_energies);     
vb_max = find(occupations==1,1,'last')

% applying scissor operator to conduction band
scissor = 0; % scissor operator in eV
band_energies_shifted = band_energies';
for i=1:nbands
    if i > vb_max
        band_energies_shifted(i) = band_energies(i) + scissor;
    end
end

% bands of interest to calculate L for 
bands_oi = [vb_max-1 vb_max vb_max+1 vb_max+2];
nbands_oi = length(bands_oi);

% calculating Lself for the bands of interst and also as function of the
% bands inlcuded in the summation
L = zeros(nbands, nbands_oi);
for count = 1:nbands_oi
    n = bands_oi(count);
    for nbands_L = 1:nbands
        for j = 1:nbands_L
            if j == n
                continue
            end
            if abs(band_energies_shifted(j)-band_energies_shifted(n)) > 0.00001
                M_z = imag(matrix_elements(n,j,1)*matrix_elements(j,n,2));      % Im{ <n|Px|j><j|Py|n> }
                L(nbands_L, count)  =  L(nbands_L, count)  -2*M_z/(band_energies_shifted(j)-band_energies_shifted(n));  
            end
            %BC(n) =  BC(n) -M_z/(band_energies(m)-band_energies(n)).^2;     % Berry curvature   
        end      
    end
end

% [L] = eV Ang^2, as it comes in VASP; we want: [L[SI]/hbar] = 1
% convert in SI units: L[1] = L[eV Ang^2] e 1E-20 m_e / hbar^2 -> same result
c = PhysicalConstants();
convert = c.e*1E-20*c.me*4*pi^2/c.h^2;
L = L*convert;

% plot the L convergence
h = plot(L,' -'); hold on; 
% set(h, {'color'}, num2cell(hsv(4),2));
xlabel('number of bands'); ylabel('L (hbar)');

% add band number as legend
for count = 1:nbands_oi
    legend_str{count} = num2str(bands_oi(count));
end
legend(legend_str)