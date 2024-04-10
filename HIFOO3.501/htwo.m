function [f,g] = htwo(A,B1,B2,C1,C2,D11,D12,D21,Ak,Bk,Ck,Dk,nvar,struct)

%   function that computes the H2 norm of the closed loop plant and its gradient
%   the plant is defined by : A,B1,B2,C1,C2,D11,D12,D21; we suppose D22 = 0
%   the controler is defined by Ak,Bk,Ck,Dk
%   the closed loop system is defined by Abig,Bbig,Cbig,Dbig 
%   f = sqrt(trace(Bbig'*X*Bbig)) where X is the solution to the Lyapunov
%   equation Abig'X + X Abig = -Cbig'Cbig
%   This routine written by Georgia Deaconu

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

n = size(A,1);
m = size(B2,2);
p = size(C2,1);

[p1,m1] = size(D11);
nk = size(Ak,1);

%matrices we are going to use to compute the gradient  
Ag = [A zeros(n,nk); zeros(nk,n) zeros(nk,nk)];
Bg = [zeros(n,nk) B2; eye(nk) zeros(nk,m)];
B1g = [B1; zeros(nk,m1)];
Cg = [zeros(nk,n) eye(nk); C2 zeros(p,nk)];
C1g = [C1 zeros(p1,nk)];
D21g = [zeros(nk,m1); D21];
D12g = [zeros(p1,nk) D12];

%the closed loop system
Abig = [A + B2*Dk*C2 B2*Ck; Bk*C2 Ak];
Bbig = [B1 + B2*Dk*D21; Bk*D21];
Cbig = [C1 + D12*Dk*C2 D12*Ck];

K = [Ak Bk; Ck Dk]; 

% %verify that the feedback system is stable
vp = eig(Abig);
absc = max(real(vp));
if (absc >= 0)
    f = inf;
    if nargout > 1
        g = nan*ones(nvar,1);
    end
    return;
end

%compute the H2 norm
[X] = lyap(Abig',Cbig'*Cbig);

f = sqrt(trace(Bbig'*X*Bbig));

if nargout > 1
    %compute the gradient
    [Y] = lyap(Abig,Bbig*Bbig');
    
    G = 2*((Bg'*X + D12g'*(D12g*K*Cg + C1g))*Y*Cg' + Bg'*X*(B1g + Bg*K*D21g)*D21g');
    G = G/(2*f);
    
    GAk = G(1:nk, 1:nk);
    GBk = G(1:nk, nk+1:nk+p);
    GCk = G(nk+1:nk+m, 1:nk);
    GDk = G(nk+1:nk+m, nk+1:nk+p);
    
    ga = GAk(:);
    gb = GBk(:);
    gc = GCk(:);
    
    g = [ga; gb; gc];
    
    %We have a reduced set of optimization variables imposed by the structure 
    if isstruct(struct)
        if (isfield(struct,'Ahat') && isfield(struct,'Bhat') && isfield(struct,'Chat'))
           
            if (~isempty(struct.Ahat))
                ga = struct.Ahat*ga;
            else
                ga = [];
            end
            
            if (~isempty(struct.Bhat))
                gb = struct.Bhat*gb;
            else
                gb = [];
            end
            if ~isempty(struct.Chat)
                gc = struct.Chat*gc;
            else
                gc = [];
            end
            
            g = [ga; gb; gc];
        end
        
        %we don't need to check for Dhat because it is included in V
        if (~isempty(struct.V))
            gtemp = GDk(:)';
            gd = gtemp*struct.V;
            g = [g; gd'];
        end
    end
end
    

