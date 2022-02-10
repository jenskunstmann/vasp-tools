function vasp_AnalyzeStructure()
% plot information about the lattice crystal and others
% use the general crystall tools routine

global SYS

file_contcar = sprintf('%s/%s', SYS.path, SYS.contcar)
cry__AnalyzeStructure(file_contcar, 'vasp');