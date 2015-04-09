function [r_i_k,r_k_alpha,I_r,l_r]=expand_field_kl_pce( r_stdnor_func, cov_r_func, pos, G_N, degPC, degKL, l_r, varargin )

options=varargin2options(varargin);
[cov_gam_func,options]=get_option( options, 'cov_gam_func', [] );
[projection_method,options]=get_option( options, 'projection_method', true );
[mean_func,options]=get_option( options, 'mean_func', [] );
[kl_eps,options]=get_option( options, 'eps', [] );
check_unsupported_options(options,mfilename);


%[r_i_alpha, I_r]=expand_field_pce_sg( r_stdnor_func, cov_r_func, cov_gam_func, pos, G_N, p_r, m_r, 'mean_func', mean_func );
%raeumlich pktweise Koeffs der PC-Entwicklung, berechnet per NORTA 
% Step 1: calculate the rho_k(pos) (Koeffs der pce, non-spatial) numerically
%im orig: [p_trans,options]=get_option( options, 'p_trans', min(pcDeg,7) );
%Warum?
rho_k=pce_expand_1d(rho_stdnor_func,degPC);
if m_gam==0
    r_j_alpha=repmat(rho_k(1), size(pos,2), 1);
    I_r=multiindex(0,0);
    C_r=ones(size(pos,2));
    return
end

% Step 0: if the number of random variables is zero

%Step 1 calculate covariance function for passing to fourier expansion
gridX = linspace(min(pos(:,1)),max(pos(:,1)),200);
gridY = linspace(min(pos(:,2)):0.02:max(pos(:,2)),200);


[X,Y] = meshgrid(gridX,gridY);
f=cos(2*pi* ((X*3)/(gridX(end)-gridX(1))+Y*3/(gridY(end)-gridY(1)))); 
[Coeff_, spatialBasis_]=expand_field_fourier2d(  f, gridX, gridY, degKL, degKL);
% ist separate Wahl v. degX u degY sinnvoll?

% Step 3: Calculate lamda_i and r_i (i.e. do Fourier (KL) expansion)
% g contains the product sqrt(lambda_i)*g_i of the KL of gamma
kl_options.correct_var=true;
%g_j_i=kl_solve_evp( C_gam, G_N, m_gam, kl_options );

% Step 5: transform gam(pos) into u
[r_j_alpha,I_r]=pce_transform_multi(spatialBasis_*Coeff, rho_k); %g_j_i, rho_k ); %trafo Koeffs rho der PC von rho_stdnor_func in KL
%aus den raeumlichen EFs ZF-vars machen

%weg:
%if projection_method
    %, 'eps' kl_eps
%    C=covariance_matrix( pos, cov_r_func );
%    v_r_i=kl_solve_evp( C, G_N, l_r );      
%    [r_i_k,r_k_alpha]=project_pce_on_kl( r_i_alpha, I_r, v_r_i );%Projektion der PC_Koeffs auf KL-Eigenfkt'en
%else
%    [r_i_k,r_k_alpha]=pce_to_kl( r_i_alpha, I_r, l_r, G_N );
%end
%r_i_k: scaled KLE eigenfunctions with mean added as first EFunction
%r_k_alpha: Coefficients w.r.t. KL EF base with some extra block due to first EV  

