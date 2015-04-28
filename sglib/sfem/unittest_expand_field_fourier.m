function unittest_expand_field_fourier
% UNITTEST_EXPAND_FIELD_FOURIER Test the EXPAND_FIELD_FOURIER and 
% EXPAND_FIELD_FOURIER2D function.
%
% Example (<a href="matlab:run_example unittest_">run</a>)
%   unittest_
%
% Calculates 1-D fourier transform
% See also EXPAND_FIELD_FOURIER, MUNIT_RUN_TESTSUITE 

%   T.Moshagen
%   Copyright 2010, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.
clear;
munit_set_function( 'expand_field_fourier' );

%[K,K_f]=compute_fft1d( );
case_indexing_fft2d();
%ftGaussianCov(1,0.25)
%ftGaussianCov(1,0.5)
%ftGaussianCov(1,1)
%ftGaussianCov(2,1)
%ftGaussianCov(3,1)
%ftGaussianCov(1,3) sigma ist faktor!
%bis jetzt kein Sinn in dieser Gliederung - kann evtl aufgegeben werden.
end

function [spatialBasis_, Coeff_]=compute_fft1d()
pos = 0:0.02:4;
%1.----------------
y = ones(size(pos,2),1); %cos(i*2*pi*pos)
abstol=15/size(pos,2)
reltol=1/size(pos,2)
if 1
deg=size(pos,2)
[ Coeff_, spatialBasis_]=expand_field_fourier(  y, pos, deg);
expected_res = zeros(size( Coeff_));
expected_res(1)=1;
assert_equals( Coeff_, expected_res, 'ft of constant','abstol', abstol,...
    'reltol', reltol);%???'fuzzy', true );
end

%2.-----------------
if 1
clear Coeff_  spatialBasis_;
tstep=0.01
%pos = tstep:tstep:1; %Periodizitaet heisst: erster u letzter Pkt sind nachbarn, nicht gleich.
%Scheinbar doch nicht: es wird mit dem abgeschlossenen Intervall gearbeite
pos = 0:tstep:1%-tstep; %aequivalent
deg=size(pos,2)-3
y=zeros(size(pos,2),1);
y = cos(2*pi*pos);
[ Coeff_, spatialBasis_]=expand_field_fourier(  y, pos, deg);

plot(pos,y),
hold on;
plot(pos, spatialBasis_*Coeff_)
expected_res = zeros(size( Coeff_));
expected_res(2)=1
assert_equals( Coeff_, expected_res, 'ft of cos(2*pi*pos)','abstol', abstol,...
    'reltol', reltol);
end

if 1
clear Coeff_  spatialBasis_;
y = cos(4*pi*pos);
[ Coeff_, spatialBasis_]=expand_field_fourier(  y, pos, deg);
plot(pos,y);
hold on;
plot(pos, spatialBasis_*Coeff_);
expected_res = zeros(size( Coeff_));
expected_res(4)=1;
assert_equals( Coeff_, expected_res, 'ft of cos(4*pi*pos)', 'abstol', abstol,...
    'reltol', reltol);
end

if 1
clear Coeff_  spatialBasis_;
y = sin(4*pi*pos);
[ Coeff_, spatialBasis_]=expand_field_fourier(  y, pos, deg);
plot(pos,y);
hold on;
plot(pos, spatialBasis_*Coeff_);
expected_res = zeros(size( Coeff_));
expected_res(5)=1;
assert_equals( Coeff_, expected_res, 'ft of sin(4*pi*pos)','abstol', abstol,...
    'reltol', reltol);
end
end



function [spatialBase_, coeff_]=case_indexing_fft2d()
gridX = 0:0.02:4;
gridY = 0:0.02:3;
degX=10
degY=10
%implemented in expand_field_fourier2d:
degX=min(degX, size(gridX,2)/2) 
degY=min(degY, size(gridY,2)/2)
%1.----------------
%xlaeuft in Zeilenrichtung, y in spaltenricht.
y = ones(size(gridY,2),size(gridX,2)); %cos(i*2*pi*pos)
abstol=15/min(size(gridX,2),size(gridY,2));
reltol=1/min(size(gridX,2),size(gridY,2));

if 0
expected_res = zeros(degY,4*degX+1);
expected_res(1,2*degX+1)=1;
[ coeff_, spatialBase_]=expand_field_fourier2d(  y,  degX, degY);

%backTrafo = sum(Coeff_.*spatialBasis_;
%surf()
assert_equals( coeff_, expected_res, '2d-ft of const f=1','abstol', abstol,...
    'reltol', reltol);%???'fuzzy', true );
assert_equals( spatialBase_(1,2*degX+1,:), ones(1,1,size(gridX,2)*size(gridY,2))...
        , '2d-ft of const','abstol', abstol...
        , 'reltol', reltol);
end

%2.------------------------------------------------------------------------
if 0
clear y Coeff_  spatialBasis_;
%centered!
gridX = -1:0.02:1;
gridY = -1:0.02:1;
[X,~] = meshgrid(gridX,gridY);
 f= sin(2*pi*(X)/(gridX(end)-gridX(1))); 
%muss laut   spatialBasis_(,,:)=reshape(cos(2*pi*(X*()+Y*(k2-1))), nPts,1); 
% Koeffs 
expected_res = zeros(degX, 4*degY+1);
expected_res(1,2*degX+1+3)=0.5;
expected_res(1,2*degX+1-3)=-0.5;
%geben:
[ coeff_, spatialBase_]=expand_field_fourier2d(  f,  degX, degY);

backTrafo=inverseFourier(coeff_, spatialBase_);
surf(reshape(backTrafo,size(X)));
assert_equals( coeff_, expected_res, '2d-ft of sin(2*pi*(X)/(gridX(end)-gridX(1))','abstol', abstol,...
    'reltol', reltol);%???'fuzzy', true );
assert_equals( spatialBase_(1,2*degX+1+3,:), reshape(f, 1,1,[])...
        , '2d-ft of sine','abstol', abstol...
        , 'reltol', reltol);
end

%3.------------------------------------------------------------------------

if 0
clear y Coeff_  spatialBasis_;
gridX = -1:0.02:1;
gridY = -1:0.02:1;
[X,~] = meshgrid(gridX,gridY);
 f= cos(pi*2*(X*3)/(gridX(end)-gridX(1))); 
%muss laut   spatialBasis_(,,:)=reshape(cos(2*pi*(X*()+Y*(k2-1))), nPts,1); 
%einen Koeff 
expected_res = zeros(degX, 4*degY+1);
expected_res(1,2*degX+1+6)=0.5;
expected_res(1,2*degX+1-6)=0.5;
%geben:
[ coeff_, spatialBase_]=expand_field_fourier2d(  f,  degX, degY);

surf(reshape(spatialBase_(1,2*degX+7,:),size(X)))
backTrafo=inverseFourier(coeff_, spatialBase_);
surf(reshape(backTrafo,size(X)));
assert_equals( coeff_, expected_res, '2d-ft of cos(2*pi*(X*3)/(gridX(end)-gridX(1))','abstol', abstol,...
    'reltol', reltol);%???'fuzzy', true );
assert_equals( spatialBase_(1,2*degX+7,:), reshape(f, 1,1,[])...
        , '2d-ft of cosine','abstol', abstol...
        , 'reltol', reltol);
end

%4.------------------------------------------------------------------------
if 0
clear y Coeff_  spatialBasis_;
gridX = -1:0.02:1;
gridY = -1:0.02:1;
[X,Y] = meshgrid(gridX,gridY);
 f= sin(2*pi*(Y)/(gridX(end)-gridX(1))); 
%muss laut   spatialBasis_(,,:)=reshape(cos(2*pi*(X*()+Y*(k2-1))), nPts,1); 
% Koeffs 
expected_res = zeros(degX, 4*degY+1);
expected_res(2,2*degX+2)=0.5;
expected_res(2,2*degX)=0.5;
%geben:
[ coeff_, spatialBase_]=expand_field_fourier2d(  f,  degX, degY);

backTrafo=inverseFourier(coeff_, spatialBase_);
surf(reshape(backTrafo,size(X)));
assert_equals( coeff_, expected_res, 'centered 2d-ft of sin(2*pi*(Y)/(gridX(end)-gridX(1))','abstol', abstol,...
    'reltol', reltol);%???'fuzzy', true );
assert_equals( spatialBase_(2,2*degX+2,:), reshape(f, 1,1,[])...
        , '2d-ft of sine','abstol', abstol...
        , 'reltol', reltol);
end

%5.------------------------------------------------------------------------
if 0
clear y Coeff_  spatialBasis_;
gridX = -1:0.02:1;
gridY = -1:0.02:1;
[X,Y] = meshgrid(gridX,gridY);
 f= cos(2*pi*(Y)/(gridX(end)-gridX(1))); 
%muss laut   spatialBasis_(,,:)=reshape(cos(2*pi*(X*()+Y*(k2-1))), nPts,1); 
% Koeffs 
expected_res = zeros(degX, 4*degY+1);
expected_res(2,2*degX+1)= 1.0;
%geben:
[ coeff_, spatialBase_]=expand_field_fourier2d(  f,  degX, degY);

backTrafo=inverseFourier(coeff_, spatialBase_);
surf(reshape(backTrafo,size(X)));
 assert_equals( coeff_, expected_res, 'centered 2d-ft of cos(2*pi*(Y)/(gridX(end)-gridX(1))','abstol', abstol,...
    'reltol', reltol);%???'fuzzy', true );
assert_equals( spatialBase_(2,2*degX+1,:), reshape(f, 1,1,[])...
        , '2d-ft of cos(2*pi*(Y)/(gridX(end)-gridX(1))','abstol', abstol...
        , 'reltol', reltol);
assert_equals( backTrafo, reshape(f, [],1)...
        , 'Backtrafo of 2d-ft of cos(2*pi*(Y)/(gridX(end)-gridX(1))','abstol', abstol...
        , 'reltol', reltol);
end

%6.-eher redundant-----------------------------------------------------------------------
if 0
clear y Coeff_  spatialBasis_;
gridX = -1:0.02:1;
gridY = -1:0.02:1;
[X,Y] = meshgrid(gridX,gridY);
 f= cos(4*pi*(Y)/(gridX(end)-gridX(1))); 
%muss laut   spatialBasis_(,,:)=reshape(cos(2*pi*(X*()+Y*(k2-1))), nPts,1); 
% Koeffs 
expected_res = zeros(degX, 4*degY+1);
expected_res(3,2*degX+1)= 1.0;
%geben:
[ coeff_, spatialBase_]=expand_field_fourier2d(  f,  degX, degY);

backTrafo=inverseFourier(coeff_, spatialBase_);
surf(reshape(backTrafo,size(X)));
 assert_equals( coeff_, expected_res, 'centered 2d-ft of cos(4*pi*(Y)/(gridX(end)-gridX(1))','abstol', abstol,...
    'reltol', reltol);%???'fuzzy', true );
assert_equals( spatialBase_(3,2*degX+1,:), reshape(f, 1,1,[])...
        , '2d-ft of cos(4*pi*(Y)/(gridX(end)-gridX(1))','abstol', abstol...
        , 'reltol', reltol);
assert_equals( backTrafo, reshape(f, [],1)...
        , 'Backtrafo of 2d-ft of cos(4*pi*(Y)/(gridX(end)-gridX(1))','abstol', abstol...
        , 'reltol', reltol);
end

if 1
clear y f Coeff_  spatialBasis_;
gridX = -1:0.02:1;
gridY = -1:0.02:1;
[X,Y] = meshgrid(gridX,gridY);
 f=cos(2*pi* ((X)/(gridX(end)-gridX(1))+Y/(gridY(end)-gridY(1)))); 
%muss laut   spatialBasis_(k1, 2*k2-1,:)=reshape(cos(2*pi*(X*(k1-1)+Y*(k2-1))), nPts,1); 
%einen Koeff 
expected_res = zeros(degX, 2*degY);
expected_res(2,2*degX+3)=0.5;
expected_res(2,2*degX-1)=0.5;
%geben:
[ coeff_, spatialBase_]=expand_field_fourier2d(  f, degX, degY);

backTrafo=inverseFourier(coeff_, spatialBase_);
surf(reshape(backTrafo,size(X)));

assert_equals( coeff_, expected_res, '2d-ft of cos(2*pi* ((X*3)/(gridX(end)-gridX(1))+Y*3/(gridY(end)-gridY(1))))'...
    ,'abstol', abstol*1.5,...
    'reltol', reltol*1.5);
assert_equals( spatialBase_(4,7,:), reshape(f, 1,1,[])...
        , '2d-ft of cosine','abstol', abstol...
        , 'reltol', reltol);
assert_equals( backTrafo, reshape(f, [],1)...
        , 'Backtrafo of 2d-ft of cos(2*pi* ((X)/(gridX(end)-gridX(1))+Y/(gridY(end)-gridY(1))))','abstol', abstol...
        , 'reltol', reltol);
end
end

function ftGaussianCov(L,sigma)
gridX = -2:0.02:1.98;
gridY = -2:0.02:1.98;
%gridX = -2:0.1:2;
%gridY = -2:0.1:2;
degX=30;
degY=30;
[X,Y] = meshgrid(gridX,gridY);
 f=gaussian_covariance([reshape(X,[],1), reshape(Y,[],1)]', zeros(2,size(X,1)*size(X,2)),L, sigma);
%muss %auch eine Gauss-Fkt ergeben
f = reshape(f,size(X));
surf(f);


%geben:
[ coeff_, spatialBase_]=expand_field_fourier2d(  f, gridX, gridY, degX, degY);
surf(coeff_)
backTrafo=inverseFourier(coeff_, spatialBase_);
backTrafo = reshape(backTrafo,size(X));
surf(backTrafo);

abstol=15/min(size(gridX,2),size(gridY,2));
reltol=1/min(size(gridX,2),size(gridY,2));

%assert_equals( coeff_, expected_res, '2d-ft of cos(2*pi* ((X*3)/(gridX(end)-gridX(1))+Y*3/(gridY(end)-gridY(1))))','abstol', abstol,...
%    'reltol', reltol);%???'fuzzy', true );
assert_equals( backTrafo, f...
        , '2d-ft of gauss kernel','abstol', abstol...
        , 'reltol', reltol);
end


function [backTrafo]=inverseFourier(coeff, spatialBase)
backTrafo=zeros(size(spatialBase,3),1);%stumpf ist trumpf
for iSpat=1:size(spatialBase,3)
    for iX=1:size(coeff,1)
        for iY=1:size(coeff,2)
            backTrafo(iSpat) = backTrafo(iSpat) + coeff(iX,iY)*spatialBase(iX,iY,iSpat);
        end
    end
end
end

function [K,K_f]=compute_operators( p_k, m_f, p_f, p_u, lex_sort )

% geometry and probability distributions will be will be the same for all
% test cases since it doesn't matter much
N=50;
[pos,els,bnd_nodes]=create_mesh_1d( 0, 1, N );
G_N=mass_matrix( pos, els );
stiffness_func={@stiffness_matrix, {pos, els}, {1,2}};

stdnor_k={@gendist_stdnor, {'beta', {2,4}, 0.001, 1.0 }};
cov_k={@exponential_covariance,{0.05,1}};

stdnor_f={@gendist_stdnor, {'beta', {2,4}, 0.001, 1.0 }};
cov_f={@exponential_covariance,{0.05,1}};


% expand the field
%expand_field_fourier_sg(  func,... %cov_gam_func, 
%pos, deg)
%[k_i_k,k_k_alpha,I_k]=expand_field_kl_pce( stdnor_k, cov_k, pos, G_N, p_k, m_k, 1 );

[I_k,I_f,I_u]=multiindex_combine({I_k,I_f},p_u);
if lex_sort
    I_u=sortrows(I_u);
end

K=kl_pce_compute_operator(k_i_k, k_k_alpha, I_k, I_u, stiffness_func, 'tensor');
K_f=kl_pce_compute_operator_fast(k_i_k, k_k_alpha, I_k, I_u, stiffness_func, 'tensor');
end