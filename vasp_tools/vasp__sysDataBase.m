function mysys = vasp__sysDataBase(ID)

% This file contains a data base of all systems you want to examine with 
% VASP TOOLS. Each data base entry contains system specific data that
% cannot be read from the VASP output files. These are filenames, K-point
% labels, the energy window for plotting, etc.
% 
% For every calculation you do (each system you study) you have to specify
% a separate entry in the list below. Each entry must contain the following
% elements:

% sys(n).ID       = some string that describes the system
% sys(n).path     = path of the directory that contains all relevant files,
%                   WITHOUT the terminal '/' (UNIX) or '\' (Windows)
% sys(n).eigenval = name of the relevant EIGENVAL file
% sys(n).ibzkpt   = IBZKPT file used to prepare self-consistent band calculation
% sys(n).ebs      = effective band structure .dat file as produced by BandUP
% sys(n).procar   = name of the relevant PROCAR file
% sys(n).doscar   = name of the relevant DOSCAR file 
% sys(n).klabels  = cell array containing the label strings of the special
%                   points; MAKE SURE THAT THE NUMBER OF LABELS IS CORRECT
% sys(n).emin     = energy range for plotting DOS and band structures in eV
% sys(n).emax     = energy range ...
% sys(n).lorbit   = VASP LORBIT parameter used to calculate the PROCAR
%                   file, defines its format

% path to directory that holds the VASP data
respath = 'F:\tud\res';

% initialization 
n=0;


%%%%%%%%%%%%%%%%%%%%%
%%%% NEW DATASET %%%%  %%%%%%%%%% GaAs %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
n=n+1;
sys(n).ID       = 'GaAs-PBE';  
sys(n).path     = [respath '/GaAs/PBE'];
sys(n).contcar  = 'static/CONTCAR';
sys(n).procar   = 'bands/PROCAR';
sys(n).eigenval = 'bands/EIGENVAL';
sys(n).doscar   = 'dos/DOSCAR';
sys(n).klabels  = {'L', '$\Gamma$','X', 'U/K', '$\Gamma$'}; 
sys(n).emin     = -10; %-2; 
sys(n).emax     = 10; %2.5; 
sys(n).lorbit   = 11;
%sys(n).cofz     = 10; % center of mass of 2L in z diredction

%%%%%%%%%%%%%%%%%%%%%
%%%% NEW DATASET %%%%  
%%%%%%%%%%%%%%%%%%%%%
n=n+1;
sys(n).ID       = 'GaAs-PBE-SOC';  
sys(n).path     = [respath '/GaAs/PBE+SOC'];
sys(n).contcar  = 'static/CONTCAR';
sys(n).procar   = 'bands/PROCAR';
sys(n).eigenval = 'bands/EIGENVAL';
sys(n).doscar   = 'dos/DOSCAR';
sys(n).klabels  = {'L', '$\Gamma$','X', 'U/K', '$\Gamma$'}; 
sys(n).emin     = -2; 
sys(n).emax     = 2.5; 
sys(n).lorbit   = 11;

%%%%%%%%%%%%%%%%%%%%%
%%%% NEW DATASET %%%%  
%%%%%%%%%%%%%%%%%%%%%
n=n+1;
sys(n).ID       = 'GaAs-HSE06';  
sys(n).path     = [respath '/GaAs/HSE06'];
sys(n).contcar  = 'static-2/CONTCAR';
%sys(n).procar   = '';
sys(n).eigenval = 'bands-1/EIGENVAL';
sys(n).ibzkpt   = 'static-2/IBZKPT';
sys(n).doscar   = 'static-7/DOSCAR';
sys(n).klabels  = {'L', '$\Gamma$','X', 'U/K', '$\Gamma$'}; 
sys(n).emin     = -2; 
sys(n).emax     = 3; 
sys(n).lorbit   = 11;

%%%%%%%%%%%%%%%%%%%%%
%%%% NEW DATASET %%%%  
%%%%%%%%%%%%%%%%%%%%%
n=n+1;
sys(n).ID       = 'GaAs-HSE06+SOC';  
sys(n).path     = [respath '/GaAs/HSE06+SOC'];
sys(n).contcar  = 'static/CONTCAR';
%sys(n).procar   = '';
sys(n).eigenval = 'bands/EIGENVAL';
sys(n).ibzkpt   = 'static/IBZKPT';
sys(n).doscar   = 'static-dos/DOSCAR';
sys(n).klabels  = {'L', '$\Gamma$','X', 'U/K', '$\Gamma$'}; 
sys(n).emin     = -5; %-2.5; 
sys(n).emax     = 5; %3.5; 
sys(n).lorbit   = 11;

%%%%%%%%%%%%%%%%%%%%%
%%%% NEW DATASET %%%%  
%%%%%%%%%%%%%%%%%%%%%
n=n+1;
sys(n).ID       = 'GaAs-HSE06+semicore';  
sys(n).path     = [respath '/GaAs/HSE06+semicore'];
sys(n).contcar  = 'static/CONTCAR';
%sys(n).procar   = '';
sys(n).eigenval = 'bands/EIGENVAL';
sys(n).ibzkpt   = 'static/IBZKPT';
sys(n).doscar   = 'static-dos/DOSCAR';
sys(n).klabels  = {'L', '$\Gamma$','X', 'U/K', '$\Gamma$'}; 
sys(n).emin     = -2; 
sys(n).emax     = 3; 
sys(n).lorbit   = 11;


%%%%%%%%%%%%%%%%%%%%%
%%%% NEW DATASET %%%%  
%%%%%%%%%%%%%%%%%%%%%
n=n+1;
sys(n).ID       = 'GaAs-HSE06+SOC+semicore';  
sys(n).path     = [respath '/GaAs/HSE06+SOC+semicore'];
sys(n).contcar  = 'static/CONTCAR';
%sys(n).procar   = '';
sys(n).eigenval = 'bands/EIGENVAL';
sys(n).ibzkpt   = 'static/IBZKPT';
sys(n).doscar   = 'static/DOSCAR';
sys(n).klabels  = {'L', '$\Gamma$','X', 'U/K', '$\Gamma$'}; 
sys(n).emin     = -2; 
sys(n).emax     = 3; 
sys(n).lorbit   = 11;

%%%% FIND DATA  %%%%
for i=1:n
    if strcmp(sys(i).ID, ID)  
        mysys = sys(i);        
        return;
    end;
end;

% print error message if data set does not exist
error(' -> SYS_ID does not exist. <-');

