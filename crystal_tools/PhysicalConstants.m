function [const] = PhysicalConstants()
%provide the elementaty physical constants in SI units

% fixed natural constantes (2018)
const.c   = 299792458; % m/s
const.h   = 6.62607015E-34; % Js = Nms = VAs^2
const.e   = 1.602176634E-19; % As = C
const.k   = 1.380649E-23; % J/K
const.NA  = 6.02214076E23; % 1/mol Avogadro constant
const.dnu = 9192631770; % 1/s frequency of hyperfine transition in 133^Cs
const.Kcd = 683; % lm/W
% derived, not yet updated, CODATA (2014), http://dx.doi.org/10.1103/RevModPhys.88.035009
const.me  = 9.10938356E-31; % kg  
const.u   = 1.660539040E-27;    % atomic mass unit [kg] 
const.eps0 =  8.854187817E-12;  % ectric constant/vacuum permittivity [As/Vm = F/m] 
const.a0  = 0.52917721067E-10;  % Bohr radius [m] 
const.Ry  = 4.359744650E-18/2;    % Rhdbger energy [J]  
%const.H  = 

end

