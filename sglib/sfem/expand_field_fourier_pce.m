function [r_i_k,r_k_alpha,I_r,l_r]=expand_field_kl_pce( rho_stdnor_func, cov_r_func, pos, G_N, degPC, degKL, l_r, varargin )

options=varargin2options(varargin);
[cov_gam_func,options]=get_option( options, 'cov_gam_func', [] );
[projection_method,options]=get_option( options, 'projection_method', true );
[mean_func,options]=get_option( options, 'mean_func', [] );
[kl_eps,options]=get_option( options, 'eps', [] );
check_unsupported_options(options,mfilename);


% check input parameters
check_condition( rho_stdnor_func, 'isfunction', false, 'rho_stdnor_func', mfilename );
check_condition( cov_r_func, 'isfunction', false, 'cov_r_func', mfilename );
check_condition( cov_gam_func, 'isfunction', true, 'cov_gam_func', mfilename );
check_range( size(pos,1), 1, 3, 'sizeof(pos,1)', mfilename );
check_condition( G_N, 'square', true, 'G_N', mfilename );
check_condition( {pos, G_N}, 'match', true, {'pos', 'G_N'}, mfilename );
check_range( degPC, 1, 10, 'degPC', mfilename );
check_range( degKL, 0, 1000, 'degKL', mfilename );
%Frage: Inwieweit macht es Sinn, die Grade G_N und degKL verschieden zu
%waehlen?
% get options
options=varargin2options( varargin );
[transform_options,options]=get_option( options, 'transform', {'correct_var', true} );
[kl_options,options]=get_option( options, 'kl_options', struct() );
%??[p_trans,options]=get_option( options, 'p_trans', min(p,7) );
[mean_func,options]=get_option( options, 'mean_func', [] );
check_unsupported_options( options, mfilename );


%[r_i_alpha, I_r]=expand_field_pce_sg( r_stdnor_func, cov_r_func, cov_gam_func, pos, G_N, p_r, m_r, 'mean_func', mean_func );
%raeumlich pktweise Koeffs der PC-Entwicklung, berechnet per NORTA 
% Step 1: calculate the rho_k(pos) (Koeffs der pce, non-spatial) numerically
%im orig: [p_trans,options]=get_option( options, 'p_trans', min(pcDeg,7) );
%Warum? degKL formerly m_gam
rho_k=pce_expand_1d(rho_stdnor_func,degPC);
if degKL==0
    r_j_alpha=repmat(rho_k(1), size(pos,2), 1);
    I_r=multiindex(0,0);
    C_r=ones(size(pos,2));
    return
end

% Step 0: if the number of random variables is zero

% abhaengig vom Eingangsgitter wird das Gitter fuer die Fouriertrafo
% erzeugt:
%Debug
%gridX = linspace(0,1,15);
%gridY = linspace(0,1,15);
%one has to take care for the domain being symmetric w.r.t. the two axes,
%otherwise cov function will not and FT will have imaginary contributions
gridX = linspace(min(pos(1,:)),max(pos(1,:)),size(pos,2));%2*size(pos,2));
gridY = linspace(min(pos(2,:)),max(pos(2,:)),size(pos,2));


% Step 2: calculate <gam_i gam_j> from <u_i u_j>
C_r=covariance_func_complete( gridX, gridY, cov_r_func );
if ~isempty( cov_gam_func )
    C_gam=covariance_matrix( posReg, cov_gam_func );
else
    if isstruct(transform_options)
        transform_options=struct2options(transform_options);
    end
    C_gam=transform_covariance_pce( C_r, rho_k, transform_options{:} );
end  
%umstrukturieren der

[EV_C_Gam, spatialBasis_]=expand_field_fourier2d(  C_gam, gridX, gridY, degKL, degKL);
% ist separate Wahl v. degX u degY sinnvoll?

% Step 3: Calculate lamda_i and r_i (i.e. do Fourier (KL) expansion)
% g contains the product sqrt(lambda_i)*g_i of the KL of gamma
kl_options.correct_var=true;
%g_j_i=kl_solve_evp( C_gam, G_N, m_gam, kl_options ); bekannt

%ziehen EV in die Basis:
g_i_j=scaleBase(sqrt( EV_C_Gam), spatialBasis_);
% Step 5: transform gam(pos) into u
[r_j_alpha,I_r]=pce_transform_multi(g_i_j, rho_k); %g_j_i, rho_k ); 
%aus den raeumlichen EFs ZF-vars machen


%if projection_method
    %, 'eps' kl_eps
    C=covariance_matrix( pos, cov_r_func );
    %v_r_i=kl_solve_evp( C, G_N, l_r ); 
    [EV_C, spatialBasis_]=expand_field_fourier2d(  C, gridX, gridY, degKL, degKL);
    %in the fourier case, spatial basis stays the same for C and C_gam,
    %recalculation is redundant
    g_i_j=scaleBase(sqrt( EV_C), sqrt(spatialBasis_));
    [r_i_k,r_k_alpha]=project_pce_on_kl( r_i_alpha, I_r, g_i_j);%v_r_i );%Projektion der PC_Koeffs auf KL-Eigenfkt'en
%else
%    [r_i_k,r_k_alpha]=pce_to_kl( r_i_alpha, I_r, l_r, G_N );
%end
%r_i_k: scaled KLE eigenfunctions with mean added as first EFunction
%r_k_alpha: Coefficients w.r.t. KL EF base with some extra block due to first EV  
%I_r Multiindex m_gam KL-Basisfkt g_i_j(x) \times p PC-Bas-Fcts = 66

%product space base Fcts
% Replace the mean ???
if ~isempty(mean_func)
    r_j_mean=funcall( mean_func, pos );
    r_j_alpha(:,1)=r_j_mean;
end

function scaledBase = scaleBase(scale_fac, spatialBasis_)
scaledBase=zeros(size(spatialBasis_));
for i=1:size(scale_fac,1)
    for j=1:size(scale_fac,2)
        scaledBase(i,j,:)=spatialBasis_(i,j,:)*scale_fac(i,j);
    end
end

function C=covariance_func_complete( gridX,gridY, covar_func )
n=size(gridX,2);
C=zeros(n,n);

for i=1:n
    for j=i:n
    C(i,j)=funcall( covar_func,[0,0]', ... repmat(pos(:,i),1,n-i+1), ...
       [gridX(i) gridY(j)]');% pos(:,i:end) );
    end
    C(i,1:i-1)=C(1:i-1,i)';
end

