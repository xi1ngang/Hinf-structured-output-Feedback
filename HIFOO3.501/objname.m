function name = objname(obj, pars)
% return a string associated with objective obj
% pars.epsilon is needed if obj is 'p'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HIFOO, A Matlab package for Fixed Order H-infinity and H2 control 
%%  Copyright (C) 2010  Marc Millstone, Michael Overton, Georgia Deaconu
%%
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch obj
    case 'h'
        name = 'H-infinity performance';
    case 'r'
        name = 'inverse of complex stability radius';
    case 's'
        name = 'spectral abscissa';
    case '+'
        name = 'stabilize only';
    case 'p'
        name = sprintf('%g-pseudospectral abscissa', pars.epsilon);
    case 't'
	    name = 'H2 performance';	
    otherwise
        name = [];
end
