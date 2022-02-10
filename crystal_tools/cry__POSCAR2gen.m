function cry__POSCAR2gen(vasp_filename)
% convert an VASP POSCAR file to a DFTB *.gen file

gen_filename = sprintf('%s.gen', vasp_filename);

[crystal] = cry__readCONTCAR(vasp_filename);
cry__WriteGEN(gen_filename, crystal);