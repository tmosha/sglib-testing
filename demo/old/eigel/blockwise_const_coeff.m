%must be such that funcall( mean_func, pos ) works.
%gives randum values with means given by coeff[1,2,3] below - so means are not random.
function [val] =  blockwise_const_coeff(pos)
% FUNCREATE 
%
% Note: 
%
% Example (<a href="matlab:run_example funcreate">run</a>)
%	blockwise_const_coeff(pos(:,1:50))
%
% See also FUNCALL, FUNCALL_FUNFUN, FUNARG, ELLIPSIS, DISP_FUNC

%   
%   Copyright 2013, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

%#ok<*AGROW>
coeff1=0.1	%left square
coeff2=0.2	%upper right square
coeff3=0.3	%lower square

for i =1:size(pos, 2)
	if pos(1,i)<0
		val(i) = coeff1;
	else
		if pos(2,i)<0
			val(i)=coeff3;
        else
			val(i)=coeff2;
		end
	end
end
	
