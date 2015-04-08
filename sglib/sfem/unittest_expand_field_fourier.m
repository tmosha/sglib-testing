function unittest_expand_field_fourier
% UNITTEST_EXPAND_FIELD_FOURIER Test the EXPAND_FIELD_FOURIER function.
%
% Example (<a href="matlab:run_example unittest_">run</a>)
%   unittest_
%
% Calculates 1-D fourier transform
% See also EXPAND_FIELD_FOURIER, MUNIT_RUN_TESTSUITE 

%   Elmar Zander
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
compute_fft2d();
%bis jetzt kein Sinn in dieser Gliederung - kann evtl aufgegeben werden.
%assert_equals( K, K_f );


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



function [spatialBasis_, Coeff_]=compute_fft2d()
gridX = 0:0.02:4;
gridY = 0:0.02:3;
degX=10
degY=10
%implemented in expand_field_fourier2d:
degX=min(degX, size(gridX,2)/2) 
degY=min(degY, size(gridY,2)/2)
%1.----------------
y = ones(size(gridX,2),size(gridY,2)); %cos(i*2*pi*pos)
abstol=15/min(size(gridX,2),size(gridY,2));
reltol=1/min(size(gridX,2),size(gridY,2));
if 1
expected_res = zeros(degX, 2*degY);
expected_res(1,1)=1;
[ Coeff_, spatialBasis_]=expand_field_fourier2d(  y, gridX, gridY, degX, degY);

%backTrafo = sum(Coeff_.*spatialBasis_;
%surf()
assert_equals( Coeff_, expected_res, '2d-ft of const','abstol', abstol,...
    'reltol', reltol);%???'fuzzy', true );
assert_equals( spatialBasis_(1,1,:), ones(1,1,size(gridX,2)*size(gridY,2))...
        , '2d-ft of const','abstol', abstol...
        , 'reltol', reltol);
end

%2.-----------------
if 1
clear y Coeff_  spatialBasis_;
[X,~] = meshgrid(gridX,gridY);
 f= cos(2*pi*(X*3)); 
%muss laut   spatialBasis_(k1, 2*k2-1,:)=reshape(cos(2*pi*(X*(k1-1)+Y*(k2-1))), nPts,1); 
%einen Koeff 
expected_res = zeros(degX, 2*degY);
expected_res(4,1)=1;
%geben:
[ Coeff_, spatialBasis_]=expand_field_fourier2d(  f, gridX, gridY, degX, degY);
%backTrafo = sum(Coeff_.*spatialBasis_;
%surf()
assert_equals( Coeff_, expected_res, '2d-ft of const','abstol', abstol,...
    'reltol', reltol);%???'fuzzy', true );
assert_equals( spatialBasis_(4,1,:), reshape(f, [],1)...
        , '2d-ft of cosine','abstol', abstol...
        , 'reltol', reltol);
end

if 1
clear y f Coeff_  spatialBasis_;
[X,Y] = meshgrid(gridX,gridY);
 f=reshape(cos(2*pi*(X*3+Y*3)), nPts,1); 
%muss laut   spatialBasis_(k1, 2*k2-1,:)=reshape(cos(2*pi*(X*(k1-1)+Y*(k2-1))), nPts,1); 
%einen Koeff 
expected_res = zeros(size( Coeff_));
expected_res(4,7)=1;
%geben:
[ Coeff_, spatialBasis_]=expand_field_fourier2d(  f, gridX, gridY, degX, degY);
y = ones(size(gridX,2),size(gridY,2)); %cos(i*2*pi*pos)

%backTrafo = sum(Coeff_.*spatialBasis_;
%surf()
expected_res = zeros(size( Coeff_));
expected_res(1,1)=1;
assert_equals( Coeff_, expected_res, '2d-ft of const','abstol', abstol,...
    'reltol', reltol);%???'fuzzy', true );
assert_equals( spatialBasis_(1,1,:), ones(1,1,size(gridX,2)*size(gridY,2))...
        , '2d-ft of cosine','abstol', abstol...
        , 'reltol', reltol);
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
