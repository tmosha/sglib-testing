function rho_N = transform_cov_to_normal_fsolve(X1, X2, rho_X, quad)
% TRANSFORM_COV_TO_NORMAL_FSOLVE Compute covariance of Gaussians using fsolve.

% Set up the function to be solved for (quad determines whether the
% function using Gaussian quadrature or using Quasi Monte Carlo is used.)
if quad
    func = @(rho_N)(transform_cov_from_normal_quad(X1, X2, rho_N)-rho_X);
else
    func = @(rho_N)(transform_cov_from_normal_qmc(X1, X2, rho_N)-rho_X);
end

% turn of annoying output
opts = optimset('Display', 'off');

% Now solve (default for rho_N0 is rho_X, because often they are not very
% far apart)
rho_N0 = rho_X; 
if exist('fsolve')
    % fsolve if not builtin (it's in the optimisation toolbox) but faster
    [rho_N, ~, flag] = fsolve(func, rho_N0, opts);
else
    % use the builtin function fzero only if fsolve is not available
    [rho_N, ~, flag] = fzero(func, rho_N0, opts);
end

% If the solver did not converge, issue a warning (better error handling
% would not be bad here)
if flag~=1
    warning('fsolve had some problems');
end
