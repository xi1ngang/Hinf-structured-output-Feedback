function [f, G] = pseudospecabsc(A,epsilon)
% call:  [f, G] = pseudospecabsc(A)
% Return epsilon-pseudospectral abscissa (max real part of epsilon-pseudospectrum) 
% of A and its gradient in complex matrix space.
% The epsilon-pseudospectrum of a matrix is the set of z for which 
% sigma_min(A-zI) <= epsilon (a positive number specified in pars.epsln).
% This is the gradient of a real function wrt the real inner product.
% This matrix has the form G = uv'/v'u, where u and v are respectively 
% normalized left and right singular vectors for the corresponding sigma_min.
% Note: <G,D> = real tr G'D = real v u'/(u'v) D = real u'Dv/u'v.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HIFOO, A Matlab package for Fixed Order H-infinity and H2 control 
%%  Copyright (C) 2010  Michael Overton
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

if any(any(isnan(A)|isinf(A)))
   f = inf;
   G = nan(size(A));
   return
end
% no advantage to calling pspa with only one output arg, unlike eig
[f, z] = pspa(A,epsilon);     % Emre Mengi's version, no mex files by default
if nargout < 2
    return
end
% z could have several components, e.g. be a conjugate pair, but we only
% care about the first - do not worry about ties or nonexistence of gradient
n = length(A);
[U,S,V] = svd(A - z(1)*eye(n));
u = U(:,n);  % left singular vector
v = V(:,n);  % right singular vector
vtu = v'*u;  % could be small, or 0 in measure zero case
G = u*v'/(vtu);  % gradient in complex matrix space