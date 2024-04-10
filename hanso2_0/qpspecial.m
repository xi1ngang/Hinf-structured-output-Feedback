function [x,d,q,info] = qpspecial(G,varargin)
% Call:
% [x,d,q,info] = qpspecial(G,varargin)
%
% Solves the QP
%
%    min   q(x) = || G*x ||_2^2 = x'*(G'*G)*x
%    s.t.  sum(x)  = 1
%              x  >= 0
%
% Inputs:
%     G            -- (M x n double) matrix G, see problem above
%     varargin{1}  -- (int) maximum number of iterates allowed
%                     If not present, maxit = max(n,10) is used
%     varargin{2}  -- (int) nmaxr: maximal number of regularizations.
%                     if G'*G is close to rank-deficient, we regularize it by
%                     adding sigma*(10^s)*eye(n), with s=0,1,...,nmaxr
%                     where sigma = 1e-12.
%     varargin{3}  -- (n x 1 double) vector x0 with initial (FEASIBLE) iterate.
%                     If not present, (or requirements on x0 not met) a
%                     feasible x0 close to the smallest column in G will
%                     be used
%
% Outputs:
%     x       -- Optimal point attaining optimal value
%     d = G*x -- Optimal point attaining optimal value
%     q       -- Optimal value found = d'*d
%     info    -- Run data:
%                info(1) =
%                   0 = everything went well, q is optimal
%                   1 = maxit reached so q might not be optimal, 
%                       but it is better than q(x0), and x is feasible
%                info(2) = #iterations used
%
%
% The problem corresponds to finding the smallest vector
% (2-norm) in the convex hull of the columns of G
%
% The method is a primal active set algorithm because the
% problem is assumed small (n <~ 30)
% and since there are only simple bounds, the EQ-QP subproblems may
% be solved efficiently by reusing factorization of the EQ-KKT matrix.
%
% Works particularly well when M >> n. I.e. it should be used when
% you need the smallest vector in the convex hull of relatively few
% points (n = columns in G) in relatively high dimensional space (dim=M).
%
% Written by Anders Skajaa
%            DTU Informatics
%            Technical University of Denmark
%            andsk@imm.dtu.dk
%
%   Send comments/bug reports to Anders Skajaa, andsk@imm.dtu.dk
%   with a subject header containing the string "hanso".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Anders Skajaa
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

n     = size(G,2);
done  = 0;
tol1  = 1e-8;
tol2  = 1e-8;
tol3  = 1e-4;
sigma = 1e-12;

if nargin > 1
    maxit = varargin{1};
else
    maxit = max(2*n,10);
end
if nargin > 2
    nmaxr = varargin{2};
else
    nmaxr = 12;
end
setx0 = 0;
if nargin > 3
    x = varargin{3};
    if size(x,1) ~= n
        setx0 = 1;
        fprintf('qpsubprob: input x0 has wrong size. using default.\n');
    end
    % This checks if initial given x0 is feasible, however it takes some
    % time -- so don't use if you're sure your initial x0 is feasible.
    % ------------------------------
    %if abs(sum(x)-1) > tol1 || min(x) < 0
    %    setx0 = 1;
    %    fprintf('qpsubprob: input x0 is not feasible. using default.\n');
    %end
    % ------------------------------
else
    setx0 = 1;
end
if setx0
    % If no x0 specified, set it to an interior point
    % (a STRICTLY feasible point) that is close to the
    % smallest column in G (in the 2-norm).
    [junk,idx] = min(sum(G.^2));
    R          = 10;
    x          = ones(n,1)/(n-1+R);
    x(idx)     = R/(n-1+R);
end

H = G'*G;

e     = ones(n,1);
otn   = 1:n;
E     = eye(n);
K     = [H,-e;-e',0]; % KKT matrix including the eq-constraint
condK = rcond(K);
u     = 0;

while condK < 1e-15 && u <= nmaxr
    u     = u + 1;
    K     = K + sigma*eye(n+1);
    sigma = 10*sigma;
    condK = rcond(K);
end

Kinv = inv(K); % this is faster than using factorizations
% because K is small. In particular when
% solving the eq-con. QP below

% W is the working set (ONLY for the ineq constraints x_i >= 0)
W     = (x < tol1); % initial working set

for k = 1:maxit

    Wlist = otn(W); % current working set for ineq-constraints

    % ======== SOLVE THE EQUALITY CON. QP HERE: ========
    h   = [-H*x;0];
    nab = sum(W);                 % number of active bounds
    D   = [E(:,W);zeros(1,nab)];  % form D
    Dt  = D';                     % Compute the transpose of D
    S   = Dt*Kinv*D;              % Compute S
    s0  = Kinv*h;                 % Using Kinv, compute s0
    lam = -S\(Dt*s0);             % S is small, so this is not slow
    ds  = Kinv*D*lam;             % Compute delta s
    p   = s0 + ds;                % the full s
    % ==================================================
    p     = p(1:n);
    pnorm = norm(p);

    if pnorm < tol2 % in this, case do not step
        [minlam,idx1] = min(lam);
        if isempty(lam)
            done = 1;
        elseif minlam >= 0
            done = 1;
        end
        if done
            x = max(x, zeros(n,1));
            x = x/sum(x);
            d = G*x;
            q = d'*d;
            info(1) = 0;
            info(2) = k;
            return;
        end
        W(Wlist(idx1)) = 0;

    else % take step
        % compute step length alpha from 16.41 N&W:
        idx2     = (~W).*(p<0);
        [m,idx3] = min( -(x.*idx2) ./ (p.*idx2) );
        alpha    = min([1,m]);

        % x_{k+1} = x_k + alpha * p:
        x = x + alpha*p;

        if alpha < 1-tol3   % if there was a blocking contraint
            W(idx3) = 1;    % add it to the working set
        end

    end

end
% if we reach this point, we never converged.
% impose exact feasibility from what we have:
x = max(x, zeros(n,1));
x = x/sum(x);
% Note that we return x and d=G*x, disregarding the
% regularization that was done
info(1) = 1;
info(2) = k;
d = G*x;
q = d'*d;
