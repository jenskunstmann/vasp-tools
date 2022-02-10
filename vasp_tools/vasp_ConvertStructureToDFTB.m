function vasp_ConvertStructureToDFTB()
% reads in the CONTCAR file and generates a CONTCAR.gen file in the same
% directory

global SYS

file_contcar = sprintf('%s/%s',SYS.path,SYS.contcar)
file_gen = sprintf('%s.gen',file_contcar);
% read in and write out
crystal = cry__readCONTCAR(file_contcar);
cry__WriteGEN(file_gen, crystal)

