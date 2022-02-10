function [kpnt_pos, eval, bandchar] = vasp__readPROCAR(filename, lorbit)
% read the projected band structure data from an PROCAR file of VASP,
% depending on the VASP 'lorbit' parameter the format of the PROCAR file is
% different
% accepted formats: lorbit = 10 (l projection), 11 (lm projection), 12 (lm+phase projection)
% lorbit = -11 (like lorbit=11 with spin-orbit coupling)
%
% reading this file usually takes very long, therefore is is only read
% once and the data is then saved in a mat-file.
%
% kpnt_pos(kpoint, component): 
%         position of the k-point in reduced coordinates,  i.e. in units of
%         the reciprocal lattice vectors 
% eval(kpoint, band) : 
%         band energies at a particular k-point 
% bandchar(kpnt,band,atom,orbital): --> no spin-polarization <--
%         the orbital=(lm) and atom projected character of a band
%
% bandchar(kpnt,band,atom,orbital,magdir): --> with SOC <--
%         magdir = 1(total, for backward compatibility), total = sqrt(mx^2 + my^2 + mz^2) (as without magnetization) 
%                  2(mx), 3(my), 4(mz), magnetization direction                  
%         if the array is addressed without explicitly specifying 'magdir'
%         then the last dimension is subsummed in the 'orbital' dimension
%         e.g. bandchar(kpnt,band,atom,orbital') 
%         then orbital' = 1...(4*9) and the first 9 elements are ususal
%         projections <- 100% backward compatible !         
%
% the orbital index in 'bandchar' is defined as:
% orbital = 1(s), 2(py), 3(pz), 4(px), 5(dxy), 6(dyz), 7(dz2), 8(dxz), 9(dx2) 
% if LORBIT = 10:
% orbital = 1(s), 2(total p), 3(0), 4(0), 5(total d), 6(0), 7(0), 8(0), 9(0) 
%
% HINT: fixed three column format problem for number of bands and number of
% atoms is fixed (if the number > 999 VASP just writes '***'): read it in
% as '%*3s'

% generate filename of the mat-file
matfile = sprintf('%s.mat',filename);
fh_mat = fopen(matfile, 'r');

% does a mat-file exist ?
if (fh_mat == -1)   % if mat-file does not exists read actual PROCAR

    % open PROCAR file
    fh = fopen(filename, 'r');
    
    % tell user to wait
    disp('Reading PROCAR file for the first time: please wait a bit.')

    % skip one header line
    skipline(fh, 1);

    % read second line
    [Q,l] = fscanf(fh, '# of k-points: %i # of bands: %i # of ions: %i \n', [3 1]);
    if l ~= 3
        error('error reading file');         
    end

    nkpnts = Q(1);   % number of k-points
    nbands = Q(2);   % number of bands (eigenvalues) per k-point
    natoms = Q(3);   % number of atoms
    
    % allocate memory
    kpnt_pos = zeros(nkpnts, 3);
    eval     = zeros(nkpnts, nbands);
    if lorbit > 0   % positive values for normal calcs
        bandchar = zeros(nkpnts, nbands, natoms, 9);    
    else            % negative values for spin-polarized or SOC calcs
        bandchar = zeros(nkpnts, nbands, natoms, 9, 4);    
    end

    % loop over all kpoints
    for kpnt = 1:nkpnts
       
        % read position of k-point
        [Q,l] = fscanf(fh, ' k-point %d :  %g %g %g weight = %g ', [5 1]);
        if l ~= 5
            error('error reading file');                 
        end   
        kpnt_pos(kpnt,:) = [Q(2) Q(3) Q(4)];

    
        % loop over all bands
        for band = 1:nbands
        
            % read band energy
            [Q,l] = fscanf(fh, ' band %*3s # energy %g # occ. %g ', [2 1]);
            if l ~= 2
                 error('error reading file');   
            end            
            eval(kpnt,band) = Q(1);            
        
            skipline(fh, 1);                   
                                                  
            % the lines following now depend on 'lorbit'
            switch(lorbit)
                case 10
                    % read band characters, 
                    [Q,l] = fscanf(fh, '%*3s %g %g %g %g', [4 natoms]);
                    if l ~= 4*natoms
                        error('error reading file');                                        
                    end  
                    bandchar(kpnt,band,:,1) = Q(1,:)';  % s character                      
                    bandchar(kpnt,band,:,2) = Q(2,:)';  % p character
                    bandchar(kpnt,band,:,5) = Q(3,:)';  % d character
                    % skip the 'tot' line and one empty line
                    skipline(fh, 2);                     
                    
                case 11
                    % read band characters, decomposed to atoms as a whole block
                    [Q,l] = fscanf(fh, '%*3s %g %g %g %g %g %g %g %g %g %g', [10 natoms]);
                    if l ~= 10*natoms
                        error('error reading file');                                        
                    end  
                    bandchar(kpnt,band,:,:) = Q(1:9,:)';                       
                    % skip the 'tot' line and one empty line
                    skipline(fh, 3);                         
%%%%%%%%%%%%%%%
                case -11 % as 11 but with SOC
                    % read band characters, decomposed to atoms as a whole block
                    [Q,l] = fscanf(fh, '%*3s %g %g %g %g %g %g %g %g %g %g', [10 4*(natoms+1)]);
                    if l ~= 10*4*(natoms+1)
                        error('error reading file');                                        
                    end  
                    % total = sqrt(mx^2 + my^2 + mz^2) (as without magnetization)  
                    bandchar(kpnt,band,:,:,1) = Q(1:9,1:natoms)';                   % total 
                    bandchar(kpnt,band,:,:,2) = Q(1:9,(natoms+2):(2*natoms+1))';    % mx
                    bandchar(kpnt,band,:,:,3) = Q(1:9,(2*natoms+3):(3*natoms+2))';  % my
                    bandchar(kpnt,band,:,:,4) = Q(1:9,(3*natoms+4):(4*natoms+3))';  % mz
                  
                    % skip the 'tot' line and one empty line
                    skipline(fh, 2);                     
%%%%%%%%%%%%%%                    
                case 12
                    % read band characters, decomposed to atoms as a whole block
                    [Q,l] = fscanf(fh, '%*3s %g %g %g %g %g %g %g %g %g %g', [10 natoms]);
                    if l ~= 10*natoms
                        error('error reading file');                                        
                    end  
                    bandchar(kpnt,band,:,:) = Q(1:9,:)';   
            
                    % skip the 'tot' and 'ion' lines, the phase factor block
                    % and one empty line after the block
                    skipline(fh, 3 + 2*natoms + 1);    
                    
                    
                otherwise
                    error('no valid lorbit option'); 
            end
        
        end % band look               
    
    end % k-point loop

    % close file
    fclose(fh);

    % save data in mat-file, v7.3 is required for very big arrays
    save(matfile,'bandchar','eval','kpnt_pos', '-v7.3');
    
else    % if mat-file does exists read it
    
    fclose(fh_mat);
    % read data from mat-file
    load(matfile,'bandchar','eval','kpnt_pos');    
    
end     % end mat-file



           

