% Setup the mesh 
scale=1; dim=2;
switch dim
    case 1
        [pos,els]=create_mesh_1d(0,1,100);
    case 2
        [pos,els]=create_mesh_2d_rect(3);
end
pos = scale*pos;
G_N=mass_matrix( pos, els );
L = size(pos,2)-1; % Set number of KL modes to the number of nodes (which is the maximum possible)

% setup the covariance function and matrix
sigma = 2;
cov_g_func = funcreate(@gaussian_covariance, funarg, funarg, scale*0.3, sigma);
C=covariance_matrix( pos, cov_g_func );

% Compute the KL for the given covariance and mass matrix
[g_i_k,sigma_k,L]=kl_solve_evp( C, G_N, L );

% Construct a PCE of a Gaussian random field around that KL
m = L;
I_g = multiindex(m, 1);
theta_k_alpha = [zeros(m,1), eye(m)];
g_i_alpha = kl_to_pce([], g_i_k, sigma_k, theta_k_alpha);

% Compute covariance matrix for that field and compare the matrices
C_KL = pce_covariance(g_i_alpha, I_g);
norm(C-C_KL)/norm(C)




%now using fourier expansion. Regular grid:
gridX = -1:0.02:1;
gridY = -1:0.02:1;
[X,Y] = meshgrid(gridX,gridY);
%(:) - reshape, [] - automatischer Nullvektor
f =funcall(cov_g_func, [X(:),Y(:)]',[])
% Compute the KL for the given covariance and mass matrix
f_M=reshape(f,size(X));
surf(f_M)
[sigma_k,g_i_k]=expandFieldFourier2dCentered(f_M,L,L);
%[ coeff_, spatialBase_]=expandFieldFourier2dCentered(  f, degX, degY);
size(sigma_k)
%kl_solve_evp( C, G_N, L );

% Construct a PCE of a Gaussian random field around that KL
m = L;
I_g = multiindex(m, 1);
theta_k_alpha = [zeros(m,1), eye(m)];
g_i_alpha = kl_to_pce([], g_i_k, sigma_k, theta_k_alpha);

% Compute covariance matrix for that field and compare the matrices
C_KL = pce_covariance(g_i_alpha, I_g);
norm(C-C_KL)/norm(C)


%D = full(sum(G_N(:)));
%1 - sum(sigma_k.^2)/(D*sigma^2)

