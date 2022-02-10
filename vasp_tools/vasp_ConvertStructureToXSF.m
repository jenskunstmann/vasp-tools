function vasp_ConvertStructureToXSF()
% reads in the CONTCAR file and generates a CONTCAR.xsf (XCrysDen)file in
% the same directory

global SYS

file_contcar = sprintf('%s/%s',SYS.path,SYS.contcar)
file_gen = sprintf('%s.xsf',file_contcar);
% read in and write out
crystal = cry__readCONTCAR(file_contcar);
cry__WriteXSF(file_gen, crystal)

