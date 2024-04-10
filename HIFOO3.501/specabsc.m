function [f, G] = specabsc(A)
% call:  [f, G] = specabsc(A)
% Return spectral abscissa (max real part of eigenvalues) of A and 
% its gradient in complex matrix space 
% In the measure zero case that the gradient is not defined 
% because of a tie for the largest real part, or because of a multiple 
% eigenvalue, the tie is broken arbitrarily.

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
if nargout < 2
   f = max(real(eig(A)));
else
   [V, Lambda] = eig(A); 
   lambda = diag(Lambda);
   [f, k] = max(real(lambda)); % spectral abscissa
   lam = lambda(k); % corresponding eigenvalue
   v = V(:,k);  % and right eigenvector
   % this "alternative" way is usually fine but it is slightly less reliable
   % left eigenvector is right eigenvector of A' for conjugate eigenvalue
   %   u_alt = (A - lam*I)'\randn(n,1); 
   %   scale = u_alt'*v;
   %   u_alt = u_alt/(scale');  % so that u_alt'*v = 1
   % the following is consistently reliable
   n = length(A);
   I = eye(n);
   e = I(k,:);
   u = e/V;   % relevant row of inverse of V ("left row eigenvector")
   % might possibly have to perturb V to make it nonsingular
   perturb = 1e-16*max(max(abs(V)));
   while isnan(u)
       % fprintf('specabsc:****** matrix of right eigenvectors is singular, perturbing\n')
       V = V + perturb*randn(size(V));
       u = e/V;
       perturb = 10*perturb;
   end
   u = u'; % by convention the conjugate transpose is usually called left eigenvector
           % (right eigenvector of A' for the conjugate eigenvalue)
           % r = norm(A'*u-lam'*u); % these are very consistent in random tests
           % r_alt = norm(A'*u_alt-lam'*u_alt);  % these vary widely
   G = u*v';  % gradient in complex matrix space
end