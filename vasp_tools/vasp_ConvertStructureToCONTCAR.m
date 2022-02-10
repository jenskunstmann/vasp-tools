function vasp_ConvertStructureToCONTCAR()
% reads in the CONTCAR file and generates a CONTCAR_new file in
% the same directory
% this can be used to convert direct into Cartesian coordinates

global SYS

file_contcar = sprintf('%s/%s',SYS.path,SYS.contcar)
file_new = sprintf('%s_new',file_contcar);
% read in and write out
crystal = cry__readCONTCAR(file_contcar);
cry__WritePOSCAR(file_new, crystal)

