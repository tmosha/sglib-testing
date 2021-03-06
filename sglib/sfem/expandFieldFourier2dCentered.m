function [ coeff_, K,phase]=expandFieldFourier2dCentered( ... 
    func, degX, degY) %  varargin )
% EXPAND_FIELD_Fourier Compute the Fourier expansion of an arbitrary real(!) Fct.
%   adapted to sglib context.
%   [coeff_,  spatialBasis_]=EXPAND_FIELD_PCE_SG( func,pos, deg) ) 
%   computes the Fourier specified by the arguments 
%   FUNC: matrix of function values evaluated e.g. at [X,Y] = meshgrid(gridX,gridY);
%   f=func(X,Y); 
% to evaluate the 
%   DEGX, DEGY: number of coefficients and  basis functions returned.
%       Note:  Independently of degX and degY, all values of FUNC are considered!
%   TODO: Sparse matrix for the fourier base plus input argument determining the 
%   limit for coefiicients whose bases should be considered.
%  by (DIFFERENTLY from MATLAB Standard function fft2() by 1/N  !)
%                    N
%    X(k) =  2* 1/N   sum  x(n)*exp(-j*2*pi*<k,x>)/N), 1 <= k <= N.
%                    n=1
%   Note: The Factor 2* is due to discrete/finite summation!
%   SpatialBasis contains the fourier basis 
%    functions  b_k=exp(-2*pi*i*<k,x>) sorted into
%    ...  cos(-x)   0           1,          0,      cos x,      sin x,      sin2x... 
%    ...  cos(-x+y) sin(y)     cos y       sin y   cos(x+y)    sin(x+y)    cos(2x +y)
%    ...cos(-2x+2y) sin(2y)    cos (2y)    sin(2y) cos(2x+2y)...
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
%sowas f. symFlag: options=varargin2options(varargin);
%[cov_gam_func,options]=get_option( options, 'cov_gam_func', [] );

%check_range( size(pos,1), 1, 3, 'sizeof(pos,1)', mfilename );
nX=size(func,2);
nY=size(func,1);

%check_range( deg, 1, N, 'deg', mfilename );

degX=min(degX,floor(nX/2));
degY=min(degY,floor(nY/2));

%symmetry? norm(func-func')
%F=fft2(func);%,NFFT);

FCentered = fftshift(fftn(ifftshift(func)));
%FTransp=FCentered'
coeff_=zeros(degY,4*degX+1);
%coeff_=zeros(degX,2*degY);
%F=fftn(ifftshift(func));


nPts = (nX*nY);
if 0
surf(imag(F+F'))
surf(real(F))
surf(real(F+F'))
end

midY = floor(size(FCentered,1)/2);
midX = floor(size(FCentered,2)/2);

coeff_(:,2*degX+1)=(real(FCentered(midY+1:1:midY+degY,midX+1)+[0; FCentered(midY:-1:midY-degY+2, midX+1)]))/nPts; %konst. Coeff
for iX=1: degX
        coeff_( 1:degY,2*degX + 2*iX)  =...
                (imag(-FCentered(midY+1:1:midY+degY,midX+ iX)  +[0;FCentered(midY:-1:midY-degY+2, midX-iX+2)]))/nPts;%???+2)))/nPts; %sinus
        coeff_( 1:degY,2*degX - 2*(iX-1))  =   ... 
                (imag(-FCentered(midY+1:1:midY+degY,midX- iX+2)+[0;FCentered(midY:-1:midY-degY+2, midX+iX)]))/nPts;%???+2)))/nPts; %sinus     
        %+coeff_( 1:degY,2*degX + 2*iX);
        coeff_( 1:degY,2*degX + 2*iX+1)  =...
                (real(FCentered(midY+1:1:midY+degY,midX+ iX+1)+[0;FCentered(midY:-1:midY-degY+2, midX-iX+1)]))/nPts; %cosinus

            coeff_( 1:degY,2*degX - 2*(iX-1)-1)  =...
                (real(FCentered(midY+1:1:midY+degY,midX- iX+1)+[0;FCentered(midY:-1:midY-degY+2, midX+iX+1)]))/nPts; 
end
%Overwite the sin2py - entries...
coeff_( 1:degY,2*degX + 2)  =...
                imag(-FCentered(midY+1:1:midY+degY,midX+ 1))/nPts;%  +[0;FCentered(midY:-1:midY-degY+2, midX-iX+2)]))/nPts;%???+2)))/nPts; %sinus
coeff_( 1:degY,2*degX    )  =                  (imag(...%-FCentered(midY+1:1:midY+degY,midX- iX+2)+
                            [0;FCentered(midY:-1:midY-degY+2, midX+1)]))/nPts;%???+2)))/nPts; %sinus     

%Remark: columns 1 and end belong to cosine().

%besser: Format
%sin_rep_k = [A_k; w_k; p_k]; %pk - Phase, hier 0 (cos) oder -pi/2 (sin)
%evaluation from -1 to 1
dim=2;
K= zeros(degY,4*degX+1,dim);
phase= zeros(degY,4*degX+1,dim);
%spatialBasis_(:,1)=1;
for k1=1:1:degY
 for k2=1:1:degX
    kX = k2-degX-1;
   K(k1, 2*k2-1,:)  = [kX, (k1-1)]; %reshape(cos(pi*(xMesh*(kX)+yMesh*(k1-1))), nPts,1); 
   K(k1, 2*k2,:)    = [(kX+1),(k1-1)]; %reshape(sin(pi*(xMesh*(kX+1)+yMesh*(k1-1))), nPts,1);
   phase(k1, 2*k2,:)    = -[ 0.5,0.0]; %lets see - cosine knows only one phase shift
   
 end
 K(k1, 2*degX+1,:) = [0, (k1-1)];
 for k2=degX+1:1:2*degX %neu: v. deg+2 an
    kX = k2-degX-1;
    K(k1, 2*k2,:)  =[kX, (k1-1)];%reshape(sin(pi*(xMesh*(kX)+yMesh*(k1-1))), nPts,1);
    K(k1, 2*k2+1,:)=[kX+1, (k1-1)]; %reshape(cos(pi*(xMesh*(kX)+yMesh*(k1-1))), nPts,1); 
    phase(k1, 2*k2,:)    = -[0.5,0.0];

%   K(k1, 2*k2-1,:)=[kX, (k1-1)]; %reshape(cos(pi*(xMesh*(kX)+yMesh*(k1-1))), nPts,1); 
%   K(k1, 2*k2,:)  =[kX, (k1-1)];%reshape(sin(pi*(xMesh*(kX)+yMesh*(k1-1))), nPts,1);
%   phase(k1, 2*k2,:)    = -[0.5,0.5];
   %M_alpha(k1,:)= [k1,k2]
 end
end
%zu debug-zwecken noch nützlich:
if 0
[xMesh,yMesh] = meshgrid(linspace(-1,1,nX),linspace(-1,1,nY));%y/(y(end)-y(1)));
spatialBasis_= zeros(degY,2*degX+1,nPts);
%spatialBasis_(:,1)=1;
for k1=1:1:degY
 for k2=1:1:degX
    kX = k2-degX-1;
   spatialBasis_(k1, 2*k2-1,:)=reshape(cos(pi*(xMesh*(kX)+yMesh*(k1-1))), nPts,1); 
   spatialBasis_(k1, 2*k2,:)  =reshape(sin(pi*(xMesh*(kX+1)+yMesh*(k1-1))), nPts,1);
   %M_alpha(k1,:)= [k1,k2]
 end
 for k2=degX+1:1:2*degX+1
    kX = k2-degX-1;
   spatialBasis_(k1, 2*k2-1,:)=reshape(cos(pi*(xMesh*(kX)+yMesh*(k1-1))), nPts,1); 
   spatialBasis_(k1, 2*k2,:)  =reshape(sin(pi*(xMesh*(kX)+yMesh*(k1-1))), nPts,1);
   %M_alpha(k1,:)= [k1,k2]
 end
end
%Remark - TODO: This delivers one column too many - last Basis should be the cosine
end