function vasp_PlotPlainBands(varargin)
% make a plain bandstructure plot for VASP, 
% k-points and distances correspond to Carestian coordinates in reciprocal space
%
% USAGE: vasp_PlotPlainBands([DeltaEF])
%
% DeltaEF = change of Fermi energy, allows to align different band plots

global SYS

% allow to chage the Fermi energy
switch(nargin)
    case(1)
        DeltaEF = varargin{1};      
        
    otherwise
        DeltaEF = 0;
end

[kpnt_pos, eval_up, eval_down] = vasp__getEigenvalues();

% read PROCAR file
% file_procar = sprintf('%s/%s',SYS.path,SYS.procar)
% [kpnt_pos, eval_up, bandchar] = vasp__readPROCAR(file_procar, SYS.lorbit);
% eval_down = [1];

size(eval_up)
size(eval_down)
size(kpnt_pos)

% get Fermi energy from DOSCAR
file_doscar = sprintf('%s/%s',SYS.path,SYS.doscar)
efermi =  vasp__getEFermi(file_doscar) + DeltaEF


% get reciprocal lattice vectors (rlv) from atomic structure from CONTCAR
% [rlv] = [1/Ang]
file_contcar = sprintf('%s/%s', SYS.path, SYS.contcar)
crystal = cry__readCONTCAR(file_contcar);
rlv  = cry__GetRLV(crystal)*2*pi;  % convert to [1/Ang]

% compile the relevant data for plotting
bands.emin = SYS.emin;        % plotting range
bands.emax = SYS.emax;                
bands.klabels = SYS.klabels;          % special point labels
bands.kpnt_pos = kpnt_pos*rlv;    % Cartesian positions of the kpoints in units of [1/Ang]
%bands.kpnt_pos = kpnt_pos;        % positions in units of the rlv
bands.bchar = 0;                  % no band characters
bands.charscal = 0;               % no fatness for the band characters = fatbands
bands.charcol = 0;                % color of the fatbands
bands.plotbands = 1;

%bands.kpnt_pos


% plot bands 
if(isempty(eval_down))     % non-spin-polarized if 'eval_down' is empty
    bands.eval = eval_up - efermi;     
    bands.linespec = '-k';             
    x = vasp__plotBands(bands);
    
    % add band numbers
    %el__plotBandNumbers(bands, x, [1 0 0], [0.0 0.0])    
    %el__getBandEnergies_MoS2(bands, x);
    %el__getBandTransitions_MoS2(bands, x);
    %el__getBandTransitions(bands, x);
    %[vbandindx, cbandindx, mbandindx, bndmax_x, bndmin_x] = el__getVCBandIndices(bands, x);
    % band masses
    %band_nr = 9;
    %nr_of_points4fit = 10;
    %[m] = el__getBandMass(x, bands, band_nr, nr_of_points4fit)
    
else                             % spin-polarized
    % spin-up bands
    bands.eval = eval_up - efermi;     
    %bands.linespec = '-r';             
    bands.linespec = '-k';             
    x = vasp__plotBands(bands);
    
    % add band numbers
    %el__plotBandNumbers(bands, x, [1 0 0], [0.0 0.0]) 
    % band masses
    band_nr = 30;
    nr_of_points4fit = 4;
    [m] = el__getBandMass(x, bands, band_nr, nr_of_points4fit)

%     % spin-down bands
%     bands.eval = eval_down - efermi;   
%     bands.linespec = '-b';             
%     vasp__plotBands(bands);
end;

title(SYS.ID,'Interpreter','none');
SetFontsInFigure(20);

% save band energies in an ASCII file, can be read in xmgrace as block data
%ematrix = [x' eval_up - efermi];
%save('bands.dat', 'ematrix', '-ASCII')



