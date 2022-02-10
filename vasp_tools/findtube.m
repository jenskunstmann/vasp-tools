function [orig, ax_vec, R] = findtube(coord, param_init)
% determine the origin 'orig', axis vector 'ax', and radius 'R' of the
% tubular atomic structure given
%
% coord(at,comp) = atomic coordinates as given in the file,
%                     in Cartesian coordinates
% param_init = [orig_x_init orig_y_init ax_phi_init ax_theta_init R_init]
%              initial guess for non-linear optimization
% that fits best to the atomic coordinates 'coord'; 
%
% origin = Cartesian coordinates of the origin of the tube
% ax_vec = Cartesian coordinates of the axial vector
% R      = Radius 
%
% units are as given by the input coordinates



% the actual fitting
[param_opt param_stdev fit_stddev CORR] = fittube(coord, param_init);


% the optimized values
orig   = param_opt(1:2); 
ax_phi   = param_opt(3); % in units of pi
ax_theta = param_opt(4); 
R        = param_opt(5);

% the standard deviation of the minimized values
orig_stdev = param_stdev(1:2);
ax_phi_stdev   = param_stdev(3); 
ax_theta_stdev = param_stdev(4); 
R_stdev = param_stdev(5);

% display the values, we show the (statistical) standard deviations
fprintf('\noptimized parameters:\n\n'); 
fprintf('  orig     = [%3.3f %3.3f] +/- [%3.3f %3.3f]\n',orig, orig_stdev); 
fprintf('  ax_phi   = %3.3f +/- %3.3f degree\n',ax_phi*180/pi, ax_phi_stdev*180/pi); 
fprintf('  ax_theta = %3.3f +/- %3.3f degree\n',ax_theta*180/pi, ax_theta_stdev*180/pi); 
fprintf('  diameter = %3.5f +- %3.5f\n\n',2*R, 2*R_stdev);
fprintf('  stddev of fit = %3.5f \n\n', fit_stddev); 

% calculate the Cartesian axis vector
ax_vec = [cos(ax_phi)*sin(ax_theta) sin(ax_phi)*sin(ax_theta) cos(ax_theta)]

% show correlation matrix
%CORR




function [param_opt param_stdev fit_stddev CORR] = fittube(coord, param_init)
% this routine fits the atomic position to the tube and returns the
% optimized parameters, standard deviations, and the correlations matrix
%
% coord(atnum, component) : Cartesian coordinates of the atoms in Angstrom
% param_init(:) : initial values of all parameters that are to be fitted
%
% param_opt(:)   : optimized values of all fitting parameters 
% param_stdev(:) : standard deviations of the fitting parameters
% fit_stddev     : general standard deviation of the whole fit
% CORR(:,:)      : correlation matrix of the fitting parameters; it tells
%                  us how independent different parameters are, a value
%                  close to 1 means they are correlated, i.e. not
%                  independent 

% number of atoms
natoms = size(coord,1);

% specify a dummy function 'f = f(param)' that is maped onto 
% 'atomdist = f(coord, param)', this is to pass the atomic coordinates
f = @(param)atomdist(coord, param);

% set fitting options
options = optimset('TolFun',1e-8,'TolX',1e-8,'MaxIter',200,'Display','iter','MaxFunEvals',10000);

% now fit
[param_opt,resnorm,residual,exitflag,output,lambda,JAC] =  lsqnonlin(f,param_init,[],[],options);

% general:
% resiudual = vector of errors; residual = dist = atomdist(coord, param_opt)
% resnorm  = dot(residual,residual) = sum(residual.^2) = sum of error squares
% stddev   = standard deviation between the data points and the fitted
%            tube = quality of the fit
fit_stddev = sqrt(resnorm/natoms);

% ci(numparamfits,2) matrix; the first column contains the upper and the
% second column lower boundaries of the 95% confidence intervall = 2 x
% standard deviation of the optimal parameters
ci = nlparci(param_opt,residual,'jacobian',JAC);
param_stdev = (ci(:,2) - param_opt')/2;

% determine the covariance matrix from the Jacobian
COVM=full(inv(JAC.'*JAC));

% the standard deviations of the parameters can be obtained obtained from
% the covariance matrix: sd_i = sqrt(COVM_ii) 
% param_stdev = sqrt(diag(COVM))

% We would like to present the user the correlation matrix
%   CORR_ij = COVM_ij/sqrt(COVM_ii*COVM_jj)
% it tells us how independent different parameters are, a value close to 1
% means they are correlated, i.e. not independent
numparamfit = length(param_init);
for i=1:numparamfit
    for j=1:numparamfit
        CORR(i,j) = COVM(i,j)/sqrt(COVM(i,i)*COVM(j,j));
    end
end





function dist = atomdist(coord, param)
% calculates the distance of each atom from the tube; 
% in the language of fitting this routine is given the model parameters and
% it returns the residual vector = distance of each atom from the tube;
% during the fit the residual vector is minimized
% 
% coord(atnum, component) : Cartesian coordinates of the atoms in Angstrom
% param(:)  = array of tube parameters
%
% dist(atnum) = distance of each atom from the tube;

% oaram is defined as [orig_x orig_y ax_phi ax_theta R];

orig = [param(1) param(2) 0];
ax_phi = param(3);
ax_theta = param(4);
ax = [cos(ax_phi)*sin(ax_theta) sin(ax_phi)*sin(ax_theta) cos(ax_theta)];
R = param(5);

% initialization
natoms = size(coord,1);
dist = zeros(natoms,1); 

for atom = 1:natoms
    atompos = coord(atom,:);
    distvec = dot(ax, atompos-orig)/dot(ax,ax)*ax - atompos + orig;
    dist(atom) = norm(distvec) - R;
end