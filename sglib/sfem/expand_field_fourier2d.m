function [ coeff_,  spatialBasis_]=expand_field_fourier2d( ... %rho_stdnor_func,...
    func,... %cov_gam_func, 
x,y, degX, degY) %G_N, p, m_gam, varargin )
% EXPAND_FIELD_Fourier Compute the Fourier expansion of an arbitrary real(!) Fct.
%   adapted to sglib context.
%   [coeff_,  spatialBasis_]=EXPAND_FIELD_PCE_SG( func,pos, deg) ) 
%   computes the Fourier specified by the arguments 
%   FUNC:Vector of function values 
%   X,Y: grid points on domain of func, assumed to be equidistant, as
%   needed by MESHGRID, used to evaluate the fourier basis 
%   b_k=exp(-2*pi*i*<k,x>)
%   DEGX, DEGY: number of coefficients and  basis functions returned.
%       Note:  Independently of degX and degY, all values of FUNC are considered!
%   TODO: Sparse matrix for the fourier base plus input argument determining the 
%   limit for coefiicients whose bases should be considered.
%  by (DIFFERENTLY from MATLAB Standard function fft2() by 1/N!)
%                    N
%    X(k) =  2* 1/N   sum  x(n)*exp(-j*2*pi*<k,x>)/N), 1 <= k <= N.
%                    n=1
%   Note: The Factor 2* is due to discrete/finite summation!
%   SpatialBasis contains the functions 
%   1,          0,      cos x,      sin x,      sin2x... 
%   cos y       sin y   cos(x+y)    sin(x+y)    cos(2x +y)
%   cos (2y)    sin(2y) cos(2x+2y)...
%   evaluated 
%   at x determined by x,y. This is necessary as e.g. KL expansion returns the same and following
%   functions need it.
%   COEFFs entries correspond to the spatialBasis.
%   The inverse is then computed by spatialBasis_:coeff_, : denouncing the
%   matrix scalar product.
%   instead of
%                     N
%       x(n) = (1/N) sum  X(k)*exp( j*2*pi*(k-1)*(n-1)/N), 1 <= n <= N.
%                    k=1
%
%
% Options
%   transform_options: {'correct_var', true}
%
% Example (<a href="matlab:run_example expand_field_pce_sg">run</a>)
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

%check_condition( func, 'isfunction', false, 'cov_r_func', mfilename );

%check_range( size(pos,1), 1, 3, 'sizeof(pos,1)', mfilename );
[xMesh,yMesh] = meshgrid(x/(x(end)-x(1)),y/(y(end)-y(1)));
nX=size(x,2);
nY=size(y,2);
%check_range( deg, 1, N, 'deg', mfilename );

degX=min(degX,nX/2) 
degY=min(degY,nY/2) 
F=fft2(func);%,NFFT);
coeff_=zeros(degX,2*degY);
% 1 0
%  cos sin spaltenweise abwechselnd 
%
nPts = (nX*nY)
for i=1:degY
 coeff_(1:degX, 2*i-1)  =2*real(F(1:degX, i))./nPts; %cosinus
 coeff_(1:degX, 2*i)    =-2*imag(F(1:degX, i))./nPts; %sinus
end
coeff_(1,1)          =  real(F(1,1))/nPts;  %const. Basis
coeff_(1,2)          =  0;  %const.-Anteil hat nichts Ungerades

spatialBasis_= zeros(degX,degY,nPts);
%spatialBasis_(:,1)=1;
for k1=1:1:degX
 for k2=1:1:degY
   spatialBasis_(k1, 2*k2-1,:)=reshape(cos(2*pi*(xMesh*(k1-1)+yMesh*(k2-1))), nPts,1); 
   spatialBasis_(k1, 2*k2,:)  =reshape(sin(2*pi*(xMesh*(k1-1)+yMesh*(k2-1))), nPts,1);
   %M_alpha(k1,:)= [k1,k2]
 end
end
spatialBasis_(1, 2,:) = zeros(nPts,1);
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
 
