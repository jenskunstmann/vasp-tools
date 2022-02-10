function vasp_OpticalSelectionRules()
% calculate the optical matrix elements for circular and linear polarized 
% light and determine the degree of circular spolarization
%
% the k-point is defined by the WAVEDERF file
% the transitions are given in the vb_list and cb_list

global SYS

% read in files and data
file = 'trans_K-/WAVEDERF';   % potentially add this file location to the SYS variable
filename = sprintf('%s/%s', SYS.path, file)
[P, ~, occ] = vasp__readWAVEDERF(filename);

% find VBM
vb_max = 36; %
find(occ==1,1,'last')

% pick transitions of interest
vb_list = [35 36]; %[vb_max-1, vb_max];
cb_list = [37 38]; %[vb_max+1, vb_max+2];

% % for relative comparison
% %p0 = 2.58 + 1.99i;     % monolayer MoSe2 VB->CB Px, for relative values [eV Ang]
% p0 = 1;                 % for absolute values
% p0_squ = p0 * conj(p0);

% polarization verctors
e_plus  = 1/sqrt(2)*[1, 1i, 0];
e_minus = 1/sqrt(2)*[1, -1i, 0];
%z       = [0, 0, 1];  

% some cheap headline
disp('   vb        cb         polariz sigma_plus sigma_minus x         y        z')
for vb = vb_list
    for cb = cb_list
        
        % momentum matrix element
        % complex valued and not uniquely defined by an arbitray phase
        p =  [P(cb, vb, 1) P(cb, vb, 2) P(cb, vb, 3)];
        
        % circular polarization
        p_plus  = e_plus * p.';  
        p_minus = e_minus * p.';       
        
        % oscillator strengths
        p_plus_squ  = p_plus * conj(p_plus);
        p_minus_squ = p_minus * conj(p_minus);
        p_x_squ     = p(1) * conj(p(1));
        p_y_squ     = p(2) * conj(p(2));
        p_z_squ     = p(3) * conj(p(3));
        
%         % relative oscillator strength
%         red = (p_plus_squ/p0_squ);      % |e+*Pcv|^2/p0     
%         blue = (p_minus_squ/p0_squ);    % |e-*Pcv|^2/p0
%         black = (p_z_squ/p0_squ);       % |z*Pcv|^2/p0
        
        % degree of circular polarization
        polariz = (p_plus_squ - p_minus_squ)/(p_plus_squ + p_minus_squ);
        
        %[vb, cb, polariz, red, blue, black]
        disp([vb, cb, polariz, p_plus_squ, p_minus_squ, p_x_squ, p_y_squ, p_z_squ])
    end
end


% Tomasz code, gives identical results
% % p0 = 2.58 + 1.99i;     % monolayer MoSe2 VB->CB Px, for relative values
% p0 = 1;                 % for absolute values
% 
% e_plus=1/sqrt(2)*[1, 1i, 0];  e_minus=1/sqrt(2)*[1, -1i, 0];  z=[0, 0, 1];  % polarization versors
% 
% count=0;
% for j = 1:length(cb)
%     for i = 1:length(vb)
%         count=count+1;
%         p =  [P(cb(j),vb(i),1), P(cb(j),vb(i),2), P(cb(j),vb(i),3)];   % p(c,v) = conj(p(v,c))
%         
%         % Yu et al 2018 2D Mater. 5 035021 
%         red = sqrt((p * e_plus.' * conj(p * e_plus.')) / (p0*conj(p0)));      % |e+*Pcv|/p0     
%         blue = sqrt((p * e_minus.' * conj(p * e_minus.')) / (p0*conj(p0)));   % |e-*Pcv|/p0
%         black = sqrt((p * z.' * conj(p * z.')) / (p0*conj(p0)));              % |z*Pcv|/p0
% 
%         % Jens' code
%         P_plus = (p(1) + 1i*p(2))/sqrt(2); P_min  = (p(1) - 1i*p(2))/sqrt(2);
%         P_plus_int = P_plus * conj(P_plus); P_min_int = P_min * conj(P_min);
%         polariz = (P_plus_int - P_min_int)/(P_plus_int + P_min_int);
%         
%         %out(count,:)=[real(vb(i)), real(cb(j)), p, real(red), real(blue), real(black), real(polariz)];
%         [real(vb(i)), real(cb(j)), real(red), real(blue), real(black), real(polariz)]
%     end
% end

