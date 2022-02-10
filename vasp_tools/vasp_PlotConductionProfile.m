function vasp_PlotConductionProfile()
% calculate the conduction profile G/G_0 of a 1D system from its band
% structure and plots it;
% G/G_0(E) = N(E) is nothing else than the number of bands N at the energy E
% the band structure needs to be present along the k-line Gamma-X ONLY
% otherwise the subroutine transmission() does not work
% for 2D and 3D structures this program produces nonsense;
% works for spin-polarized calculations
% [G/G_0(E)] = total number of channels, includeing spin degeneracy (factor
% 2) for non-spin-polarized calcs
% [G/G_0(E)] = number of spin-up and spin-down channels plotted separately
% for spin-polarized calcs

global SYS

% number of energy points in the energy interval between emin and emax
Npoints = 1000;

% read EIGENVAL file
file_eigenval = sprintf('%s/%s',SYS.path,SYS.eigenval);
[kpnt_pos eval_up eval_down] = vasp__readEIGENVAL(file_eigenval);

% get Fermi energy from DOSCAR
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar);
efermi = vasp__getEFermi(file_doscar);

% plot bands 
if(size(eval_down) == [1 1])     % non-spin-polarized if 'eval_down' is empty   
    
    % calculate the transmission with the subroutine from below
    [E,Trans] = transmission(eval_up-efermi, SYS.emin, SYS.emax, Npoints);
    
    % include spin degeneracy, e.g. the number of channels is doubled
    Trans = Trans*2;

    % define plotting range for the plot
    plotmax = 1.1*max(Trans);

    % plot
    hold on
    plot(E,Trans)

else                             % spin-polarized
    % calculate the transmission with the subroutine from below
    [E,Trans_up]   = transmission(eval_up-efermi, SYS.emin, SYS.emax, Npoints);
    [E,Trans_down] = transmission(eval_down-efermi, SYS.emin, SYS.emax, Npoints);
    
    % define plotting range for the plot
    plotmax = 1.1*max(Trans_up);

    % plot
    hold on
    plot(E,Trans_up,   'b');
    plot(E,Trans_down, '--r');
    
    legend('spin up','spin down')
end;

% set plot attributes
plot([0 0], [0 plotmax], '--k');  % plot vertical line for Fermi level
title(SYS.ID,'Interpreter','none');
xlabel('Energy (eV)')
ylabel('G/G_0');
ylim([0 plotmax]);  % set plotting range    




function [E,Trans]=transmission(En,Emin,Emax,Npoints)
% the program 'Transmission' calculates the transmission
% through 1D conductor which has periodic structure (e.g., nanotube).
% En = energy bands and k-points ('X' colons and 'Kpoints' rows)
% calculated on the 1D grid of Kpoints k-points
% and calculates the transmission for the energy interval
% from Emin to Emax on the equidistant grig of Npoints.
%
% (c) Viktor Bezugly, TU Dresden, 2011.

% determine the number of bands
N=size(En,2);
Kpoints=size(En,1);

EI=Emax-Emin; % energy interval
if EI <= 0
    error('incorectly defined Emin and Emax');         
end;

% allocate the array
Trans=zeros(1,Npoints);
E = zeros(1,Npoints);

% main loop over all Npoints
for i=1:Npoints
    E(i)=Emin+i*EI/Npoints; % current energy value
    
    % loop over all bands
    for j=1:N %N
        Ncross=0;
        % loop over k-points
        for k=1:Kpoints-1
            if ( (En(k+1,j)>E(i) && En(k,j)<E(i)) || (En(k+1,j)<E(i) &&  En(k,j)>E(i)) )
                Ncross=Ncross+1;
            end;
        end;
        %Trans(i)=Trans(i)+2*Ncross; % two channels per band (spin up and down)
        Trans(i)=Trans(i)+Ncross;  % one channels per band
    end;
end

