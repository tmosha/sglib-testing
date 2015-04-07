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

[K,K_f]=compute_fft1d( );
assert_equals( K, K_f );


function [spatialBasis_, Coeff_]=compute_fft1d()
pos = 0:0.02:4;
%1.----------------
y = ones(size(pos,2),1); %cos(i*2*pi*pos)
if 1
deg=size(pos,2)
[ Coeff_, spatialBasis_]=expand_field_fourier(  y, pos, deg)
expected_res = zeros(size( Coeff_));
expected_res(1)=1;
assert_equals( Coeff_, expected_res, 'ft of constant','reltol', 1/size(pos,2));%???'fuzzy', true );
end
%2.-----------------
clear;
tstep=0.01
pos = tstep:tstep:1; %Periodizitaet heisst: erster u letzter Pkt sind nachbarn, nicht gleich.
%pos = 0:tstep:1-tstep; %aequivalent
deg=size(pos,2)-3
if 1
y=zeros(size(pos,2),1);
y = cos(2*pi*pos);
[ Coeff_, spatialBasis_]=expand_field_fourier(  y, pos, deg)

plot(pos,y),
hold on;
plot(pos, spatialBasis_*Coeff_)
end

if 1
y=zeros(size(pos,2),1);
y = cos(4*pi*pos);
[ Coeff_, spatialBasis_]=expand_field_fourier(  y, pos, deg)

plot(pos,y);
hold on;
plot(pos, spatialBasis_*Coeff_);
end
expected_res = zeros(size( Coeff_));
expected_res(2)=1
assert_equals( Coeff_, expected_res );
assert_equals( Coeff_, expected_res, 'ft of constant','reltol', 1/N);

for i=1:1
    y = cos(i*2*pi*pos)
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
