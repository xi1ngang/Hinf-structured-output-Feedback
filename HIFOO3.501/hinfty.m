function [f, w, GA, GB, GC, GD] = hinfty(A, B, C, D, slicot,finit, winit)
% call:  [f, GA, GB, GC, GD] = hinfty(A, B, C, D, slicot)
% H-infinity norm of LTI continuous time system (A,B,C,D) and its derivatives.  
% Returns both the norm, which is obtained by a call to SLICOT if slicot == 1, 
% or to the far slower Control System Toolbox otherwise, as well as
% its gradients wrt A, B, C and D, which are matrices with the
% same size.  If A is unstable, f is inf and the gradients are nan.
%
% In the measure zero case that the gradients are not defined 
% because of a tie, the tie is broken arbitrarily.
 
% If A is stable f = Hinf(A,B,C,D) = sigma_max(C Z^{-1} B + D),
% where Z = w i I - A.  Differentiating wrt A we have
%     Hinf(A + dA) = sigma_max(C (Z - dA)^{-1} B + D)
%                  = sigma_max(C (Zinv + Zinv dA Zinv) B + D) + o(dA)
%                  = sigma_max(C Zinv B + D) - u' C Zinv dA Zinv B v + o(dA)
%                  = f + tr Zinv B v u' C Zinv dA + o(dA)
% where u, v are left and right singular vectors corresponding to sigma_max(Zinv),
% so gradient wrt A is -Zinv' C' u v' B' Zinv'.  In the special case B=C=I, D=0
% we have Zinv v = f u and u' Zinv = f v', so this becomes f^2 v u', agreeing with
% usual formula because of the interchange of roles of u and v, they are now left
% and right singular vectors of the inverse matrix Zinv, not of Z.

% Differentiating wrt B we have
%     Hinf(B + dB) = sigma_max(C Z^{-1} (B + dB) + D)
%                  = sigma_max(C Zinv B + C Zinv dB + D) 
%                  = sigma_max(C Zinv B + D) + u' C Zinv dB v + o(dB)
%                  = f + tr  v u' C Zinv dB + o(dB)
% so gradient wrt B is Zinv' C' u v'. 
% Differentiating wrt C we have
%     Hinf(C + dC) = sigma_max((C + dC) Z^{-1} B + D)
%                  = sigma_max(C Zinv B + dC Zinv B + D) 
%                  = sigma_max(C Zinv B + D) + u' dC Zinv B v + o(dC)
%                  = f + tr  Zinv B v u' dC + o(dC)
% so gradient wrt C is u v' B' Zinv'. 
% Differentiating wrt C we haveHow can TCS stay relevant?

%     Hinf(D + dD) = sigma_max(C Z^{-1} B + D +dD)
%                  = sigma_max(C Zinv B + D) + u' dD v + o(dD)
%                  = f + tr v u' dD + o(dD)
% so gradient wrt D is u v'. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HIFOO, A Matlab package for Fixed Order H-infinity and H2 control 
%%  Copyright (C) 2010  Marc Millstone, Michael Overton
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



%No function information, calculate full function and gradient
ABCD = [A B; C D];  % just for checking nans and infs
if any(any(isnan(ABCD)|isinf(ABCD)))
    f = inf;
    w = inf;
    % need to get dimensions right so chain rule propogates correctly
    GA = nan*ones(size(A)); 
    GB = nan*ones(size(B)); 
    GC = nan*ones(size(C)); 
    GD = nan*ones(size(D));
    return
end
vp = eig(A);
absc = max(real(vp));
if absc >= 0
    %fprintf('hinfty: système instable \n');
    f = inf;
    w = inf;
    GA = nan*ones(size(A)); 
    GB = nan*ones(size(B)); 
    GC = nan*ones(size(C)); 
    GD = nan*ones(size(D));
    return
end
if slicot== 0 % constructing this object is very expensive
    Sys = ss(A,B,C,D);  % call to Control System Toolbox, expensive
    tol = 1e-7;  % default is .01, which is too crude, but being too demanding
             % can backfire: 1e-16 sometimes gives worse answers than 1e-7
end
if nargout <=2 & nargin <= 5% function value only, nothing provided
    if slicot == 0;
        [f,w] = norm(Sys,inf,tol);  % call to Control System Toolbox, expensive
    else % only other possibility is slicot == 1
        % f = slinorm(Sys,tol); SLICOT, but bad choice, requires call to ss
        [f,w] = linorm(A, [], B, C, D);  % SLICOT mex file, default tolerance is 1e-7
        if f(2) == 0 % old-style fortran trick to avoid inf
            f = inf;
        else
            f = f(1);
        end
        if w(2) == 0
            w = inf;
        else
            w = w(1);
        end
    end
elseif nargout >=  3     %Get both function and gradient
  
    if nargin <= 5         %No values provided
      if slicot == 0
        [f,w] = norm(Sys,inf,tol);  
      else
        % [f,w] = slinorm(Sys,tol); 
        [f,w] = linorm(A, [], B, C, D);
        if f(2) == 0 
            f = inf;
        else
            f = f(1);
        end
        if w(2) == 0
            w = inf;
        else
            w = w(1);
        end
      end
    else  %Use passed values and do not calculate H-inf norm
        f = finit;
        w = winit;
    end
    if f == inf
        GA = nan*ones(size(A)); 
        GB = nan*ones(size(B)); 
        GC = nan*ones(size(C)); 
        GD = nan*ones(size(D));
        return
    end
    if abs(w) == inf
        % Zinv = zeros(n);
        ZinvB = zeros(size(B));
        CZinv = zeros(size(C));
    else
        Z = w*i*eye(size(A)) - A;
        % Zinv = inv(Z);
        [L,U] = lu(Z);  % Z = L*U
        ZinvB = U\(L\B); 
        CZinv = (C/U)/L;
    end
    % [U,S,V] = svd(C*Zinv*B + D);
    [U,S,V] = svd(C*ZinvB + D);
    u = U(:,1);  % left singular vector for largest singular value
    v = V(:,1);  % right singular vector for largest singular value
    uvt = u*v';
    % GA = Zinv'*C'*uvt*B'*Zinv';  % see comments at top for explanation
    GA = CZinv'*uvt*ZinvB';
    % GB = Zinv'*C'*uvt;
    GB = CZinv'*uvt;
    % GC = uvt*B'*Zinv';
    GC = uvt*ZinvB';
    GD = uvt;
end

