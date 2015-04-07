function [ coeff_,  spatialBasis_]=expand_field_fourier( ... %rho_stdnor_func,...
    func,... %cov_gam_func, 
pos, deg) %G_N, p, m_gam, varargin )
% EXPAND_FIELD_Fourier Compute the Fourier expansion of an arbitrary Fct.
%   [coeff_,  spatialBasis_]=EXPAND_FIELD_PCE_SG( func,pos, deg) ) 
%   computes the Fourier specified by the arguments 
%   FUNC:Vector of function values 
%   POS: grid points on domain of func, assumed to be equidistant
%   DEG: number of coefficients and  basis functions returned.
%       Note:  Independently of deg, all values are considered!
%  by (DIFFERENTLY from MATLAB Standard function fft() by 1/N!)
%                    N
%    X(k) =  2* 1/N   sum  x(n)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.
%                    n=1
%   Note: The Factor 2* is theoretically unjustified here!
%   SpatialBasis contains the functions 1, sin x, cos x, sin2x... evaluated 
%   at x determined by pos. This is necessary as KL expansion returns the same and following
%   functions need it.
%   The inverse is then computed by
%                     N
%       x(n) = (1/N) sum  X(k)*exp( j*2*pi*(k-1)*(n-1)/N), 1 <= n <= N.
%                    k=1

%
% Options
%   transform_options: {'correct_var', true}
%
% Example (<a href="matlab:run_example expand_field_pce_sg">run</a>)
%   N=51;
%   p_k=4;
%   m_k=4;
%   stdnor_k={@lognormal_stdnor,{0.5,1}};
%   stdnor_k={@beta_stdnor,{4,2}};
%   [pos,els,bnd]=create_mesh_1d( 0, 1, N );
%   G_N=mass_matrix( pos, els );
%   for i=1:4
%     lc_k=0.5^(i-1);
%     fprintf('conv. length: %g\n', lc_k);
%     cov_k={@gaussian_covariance,{lc_k,1}};
%     [k_i_alpha, I_k]=expand_field_pce_sg( stdnor_k, cov_k, [], pos, G_N, p_k, m_k );
%     subplot(2,2,i); plot_pce_realizations_1d( pos, k_i_alpha, I_k );
%   end
% 
% See also 

%   
%   Copyright 2006, Institute of Scientific Computing, TU Braunschweig.
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

% check input parameters


%check_condition( rho_stdnor_func, 'isfunction', false, 'rho_stdnor_func', mfilename );
%check_condition( func, 'isfunction', false, 'cov_r_func', mfilename );

check_range( size(pos,1), 1, 3, 'sizeof(pos,1)', mfilename );
% nicht noetig: check_range( deg, 1, 10, 'deg', mfilename );
N=size(pos,2);
%NFFT = 2^nextpow2(N);%deg soll zweierpotenz sein - bringt nix f.
%Genauigkeit
%Ferner:einen kleineren Grad der Entwicklung an fft zu übergeben,
%vermindert die Genauigkeit heftigst.
%check_range( m_gam, 0, 1000, 'm_gam', mfilename );
deg=min(deg,N);
F=fft(func);%,NFFT);
%Intervall-Trafo:
b=pos(end)-pos(1)
coeff_=zeros(2*deg-1,1);
 coeff_(1)          =real(F(1))/N;  %const. Basis
 coeff_(2:2:end-1)  =real(F(2:deg))./N; %cosinus
 coeff_(3:2:end)    =-imag(F(2:deg))./N; %sinus
spatialBasis_= zeros(size(pos,2),size(coeff_,1));
spatialBasis_(:,1)=1;
for i=2:2:size(coeff_,1)
   spatialBasis_(:,i)=cos(2*pi*pos*(i/2)); 
   spatialBasis_(:,i+1)=sin(2*pi*pos*(i/2)); 
end
% Step 2: calculate <gam_i gam_j> from <u_i u_j>
%C_r=covariance_matrix( pos, cov_r_func );
%if ~isempty( cov_gam_func )
%    C_gam=covariance_matrix( pos, cov_gam_func );
%else
%    if isstruct(transform_options)
%        transform_options=struct2options(transform_options);
%    end
%    C_gam=transform_covariance_pce( C_r, rho_k, transform_options{:} );
%end 

% Step 3: Calculate lamda_i and r_i (i.e. do KL expansion)
% g contains the product sqrt(lambda_i)*g_i of the KL of gamma
 
