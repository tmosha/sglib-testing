function [u, solve_info, state] = undamped_spring_solve(state, p, varargin)
% ELECTRICAL_NETWORK_SOLVE Short description of electrical_network_solve.
%   ELECTRICAL_NETWORK_SOLVE Long description of electrical_network_solve.
%
% Example (<a href="matlab:run_example electrical_network_solve">run</a>)
%
% See also

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

options = varargin2options(varargin);
check_unsupported_options(options, mfilename);

assert(numel(p)==2);

m = p(1)+1.5;
k = p(2)+1.5;

u0 = state.u0;
T = state.T;
solve_info = struct();
d = state.d;
switch state.solver
    case 'direct'
        u=undamped_spring_solve_direct(u0, k, m, d, T);
    case 'numerical'
        u=undamped_spring_solve_ode(u0, k, m, d, T);
    otherwise
        error('foo:bar', 'Unknown solver: %s', solver);
end

function u=undamped_spring_solve_ode(u0, k, m, d, T)
[t,u] = ode45(@undamped_spring_ode, [0, T], u0, [], k, m, d);
assert(t(end)==T)
u = u(end,:)';

function u=undamped_spring_solve_direct(u0, k, m, d, T)
x0 = u0(1);
v0 = u0(2);
assert(k>0 && m>0 && 0*d>=0);
alpha = d/m;
D = k/m-(d/m)^2;
if true && D>0
    omega = sqrt(D);
    xt = exp(-alpha*T) * (x0 * cos(omega*T) + (v0 + alpha*x0) / omega * sin(omega*T));
    vt = exp(-alpha*T) * (v0 * cos(omega*T) - ((alpha*v0 + alpha^2*x0) / omega + x0 * omega) * sin(omega*T));
elseif D<0
    omega = sqrt(-D);
    xt = exp(-alpha*T) * (x0 * cosh(omega*T) + (v0 + alpha*x0) / omega * sinh(omega*T));
    vt = exp(-alpha*T) * (v0 * cosh(omega*T) - ((alpha*v0 + alpha^2*x0) / omega - x0 * omega) * sinh(omega*T));
else
    error('not covered yet');
end

if xt>1 && d>=0
    keyboard
end
u = [xt; vt];
% xt = x0 * cos(omega*T) + v0 / omega * sin(omega*T);
% vt = v0 * cos(omega*T) - x0 * omega * sin(omega*T);
%log10(abs(u))'

%wieso undamped?
function [dudt] = undamped_spring_ode(t, u, k, m, d)
x = u(1);
v = u(2);
a = -(2*d*v+k*x)/m;
dudt = [v; a];

