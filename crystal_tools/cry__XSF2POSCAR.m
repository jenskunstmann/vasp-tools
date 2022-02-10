function cry__XSF2POSCAR(xsf_filename)
% convert an xsf file to a VASP POSCAR file

poscar_filename = sprintf('%s_POSCAR', xsf_filename);

[crystal] = cry__readXSF(xsf_filename);
cry__WritePOSCAR(poscar_filename, crystal);