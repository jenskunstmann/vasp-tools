% initialize the 'vasp tools' programs; 
% here the system under consideration is defined as SYS_ID
% the global SYS struct is read from the sysDataBase

clear SYS_ID SYS

global SYS_ID SYS

% GaAs
%SYS_ID = 'GaAs-PBE'
%SYS_ID = 'GaAs-PBE-SOC'
%SYS_ID = 'GaAs-HSE06'
%SYS_ID = 'GaAs-HSE06+SOC';
%SYS_ID = 'GaAs-HSE06+semicore'; 
SYS_ID = 'GaAs-HSE06+SOC+semicore'

SYS = vasp__sysDataBase(SYS_ID);