function K=pdetool_stiffness_matrix( pos, els, k, varargin )
% PDETOOL_STIFFNESS_MATRIX Compute stiffness matrix using the PDE Toolbox.
%   K=PDETOOL_STIFFNESS_MATRIX(POS, ELS, K) computes the stiffness matrix
%   from the mesh specified in POS and ELS using the PDE Toolbox and the
%   coefficient field specified in K.
%   Options:
%   coeff_el: true - Coefficient field evaluated at simplex centers
%               false - coefficient field given at grid points
%
% Notes:
%   This function is only available if the PDE Toolbox is available and on
%   the matlab path.
%
% Example (<a href="matlab:run_example pdetool_stiffness_matrix">run</a>)
%   [pos,els]=load_pdetool_geom('circle', 'numrefine', 0);
%   k = 3+sum(pos.^2,1)';
%   K=pdetool_stiffness_matrix(pos, els, k);
%   spy(K)
%
%   K2=stiffness_matrix( pos, els, k );
%   norm( K-K2, 'fro' )
%
% See also PDETOOL_MASS_MATRIX, STIFFNESS_MATRIX

%   Elmar Zander
%   Copyright 2013, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

if nargin<3
    error('sglib:pdetool_stiffness_matrix', 'Not enough arguments. You''re probably using the old version of the function interface');
end

options=varargin2options(varargin);
[coeff_el,options]=get_option( options, 'coeff_el', 0 );
check_unsupported_options(options, mfilename);

% pde toolbox needs mean values in the center of each element
t=[els; zeros(1, size(els,2))];
if(coeff_el)
    if ~(length(k)== size(els, 2))
        %exception
    end
    ct=k';
else
     if ~(length(k)== size(pos, 2))
     end
     cp=k';
    ct=sum(cp(t(1:3,:)),1)/3;
end
       
p=pos;


% call pde toolbox function
K=assema(p,t,ct,0,0);
