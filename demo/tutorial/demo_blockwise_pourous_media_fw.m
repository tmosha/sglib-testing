
function demo_blockwise_pourous_media_fw()
clear

% load or create the geomatry
% !!! pos, els, G_N (kann auch stattdessen I sein) und stiffness_func
% werden von dem FEM code benoetigt, stiffness_func muss so sein, dass
pos = [-1 -1
   -0.5  -1
    0   -1
    0.5 -1
    1   -1
    -1 -0.5
   -0.5  -0.5
    0   -0.5
    0.5 -0.5
    1   -0.5
   -1   0
   -0.5 0
    0   0
    0.5 0
    1   0
   -1   0.5
   -0.5 0.5
    0   0.5
 %   0.5 0.5
 %   1   0.5
    -1  1
   -0.5  1
    0   1
%needed for inflow of concentration:    
     -0.75 0.75
     -0.75 1
 %   0.5 1
%    1   1
    ]'
bval0=1.0;
bval1=0.0;
bnd_type_val = [0 0
   -0  0
    0  0
    0 0
    1   bval1
    0 0
   0  0
    0   0
    0   0
    1   bval1
    0   0
    0   0
    0   0
    0   0
    1   bval1
    0   0
    0   0
    0   0
 %   0.5 0.5
 %   1   0.5
    1   bval0
    1   bval0
    1   bval0
%needed for inflow of concentration:    
    0   0
    1   bval0
 %   0.5 1
%    1   1
    ]
els = [1 6 7 3
    1 7 2 3
    2 7 8 3
    2 8 3 3
    3 8 9 2
    3 9 4 2
    4 9 10 2
    4 10 5 2
    6 11 12 3
    6 12 7  3
    7 12 13  4
    7 13 8   4
    8 13 14  4
    8 14 9   4
    9 14 15  4
    9 15 10  4
    11 16 17 1
    11 17 12 1
    12 17 18 1
    12 18 13 1
% instead of
%   16 19 20 2
%   16 20 17 2
    %local refinement for inflow condition on concentration:
    16 19 22 1
    19 23 22 1
    23 20 22 1
    20 17 22 1
    17 16 22 1
    %
    17 20 21 1
    17 21 18 1
    ]'
%Koord-trafo auf X_min,X_max - todo! 
%eben unnoetig: G_N = pdetool_mass_matrix(pos, els);
%stiffness_func=funcreate(@pdetool_stiffness_matrix, pos, els, @funarg);
%[d,N]=size(pos);
%calculate realization of field K with const coeffs
%Klon v. test_MCMC
% UQ setup
variance_Y        = 4/32;         % variance of Y
mean_Y            = -10;         % mean value of Y

% BU setup: synthetic truth
zonal_Y_true = mean_Y + [...
%war fuer die ausgesparte Ecke:  -15;...
  - 0;...
  - 2;...
  - 1;...
  + 3];
n_zones = length(zonal_Y_true)

% define measurement locations in domain coordinates
measloc_hx = [];
measloc_hy = [];
measloc_cx = 9.0;
measloc_cy = 2.5;
measloc_Qx = [];
measloc_Qy = [];

test_flag = 0 
if test_flag
    coeff= ones(size(els,2),1);
else
    %war mal Klon v. Wolfgang
      zonal_Y_real          = zonal_Y_true + randn(n_zones,1)*sqrt(variance_Y)/10;
      coeff                = MC_write_zonal_Y(els,zonal_Y_real)%,mean_Y,  Z, n_zones);

end
%solving for head:
%assemble stiffness matrix
Mat=pdetool_stiffness_matrix( pos, els, coeff, 'coeff_el', true )
%ass. rhs ( dirichlet pressure treated below):
f=zeros(size(pos,2),1);
%Dirichlet-values:
bnd_h = find(bnd_type_val(:,1));
[P_I,P_B]=boundary_projectors( bnd_h, size(pos,2) )
% Ks=I_I*K*I_I+I_B;         % modify operator
%   fs=I_I*(f-K*I_B*g)+I_B*g; % modify RHS
%   u=Ks\fs;                  % now solve with mod operator
[Kn, fn] = apply_boundary_conditions_system( Mat, f, bnd_type_val(:,2), P_I, P_B )
%solve lin.sys

h=Kn\fn;

plot_field(pos, els(1:3,:), h);%, opts{:} ); 
%Woher kennt h die Dirichlet-Werte? Reassembliert
%evaluate measurement

end

%von Wolfgang, aber geaendert
function [Y] = MC_write_zonal_Y(els,zonal_Y)%,mean_Y,Z,n_zones)
% MC_WRITE_ZONAL_Y writes zone-wise constant vales from ZONAL_Y into
% a Y-field. Z contains zonation information for N_ZONES zones.
% NEL is number of elements in output.
% L-Shape obstacle is buit in, set to MEAN_Y-10
% 
% version 20.05.2014 / WN

%scheinbar default: Y                 = ones(nel,1)*mean_Y;
for i=1:size(els,2)
  Y(i,1)         = zonal_Y(els(4, i));
end
%for i=1:n_zones
%  Y(Z{i})         = zonal_Y(i);
%end
% enforing that the L-shape obstacle is there!
%zonal_Y(1)        = mean_Y-10;
%Y(Z{1})           = zonal_Y(1);
end
