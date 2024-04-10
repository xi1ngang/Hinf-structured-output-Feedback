function g = chainABCD(n, nhat, B2, C2, D12, D21,D22, GA, GB, GC, GD,x,structure)
% EDITED BY MARC MILLSTONE FOR D22 !=0
% Construct the gradient with respect to the parameterization specified
% in function getABCDbig, via chain rule.
% The complex gradients of the objective function of Abig, Bbig, Cbig, Dbig
% are given by GA, GB, GC, GD respectively. In principle, the objective
% function could be any function, but we are interested in two cases: the
% H-infinity norm (see hinfty.m) and the spectral abscissa (see specabsc.m;
% in this case GB, GC and GD are all trivially zero).

% Derivation of chain rule is from first principles: 
% f(Q + M(X+DeltaX) N) is approx f(Q  + MXN) + <del f, M DeltaX N> and the 
% latter is real tr G' M DeltaX N = real tr N G' M DeltaX = <M' G N', DeltaX>,
% where G is the gradient of f in matrix space.  Thus the gradient wrt X
% of f(Q + MXN) is M'GN'; equivalently, the gradient wrt vec(X) is vec(M'GN').

%for hinfty norm
%f(A,B,C,D) where each depends on matrices Ahat, Bhat, Chat, Dhat
% This gives rise to a nonlinear dependence to the above chain
% rule is a little more complicated.  First one must determin
%how A(Ahat+delta, Bhat, Chat, Dhat) = A + deltaA
%etc.  


% Note: the reason we take the real part is not because of the real inner
% product on complex matrix space.  It's because the parameter space is always
% real, even though the matrix functions may be complex.  Usually the matrix
% functions are also real, but of course the eigenvalues are generally complex.

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

indx1 = 1:n;
indx2 = n + 1 : n + nhat;

%If D22 == 0, use the old calculations, they are faster

if(max(max(abs(D22))) ==0)
    delAbigAhat = real(GA(indx2,indx2));         % [0  I]*GA*[I; 0]
    delAbigBhat = real(GA(indx2,indx1)*C2');     % etc
    delAbigChat = real(B2'*GA(indx1,indx2));     % etc
    delAbigDhat = real(B2'*GA(indx1,indx1)*C2'); % [B2' 0]*GA*[C2'; 0]
    delBbigBhat = real(GB(indx2,:)*D21');     % GB*[0; D21']     (delBbigAhat is 0)                                          
    delBbigDhat = real(B2'*GB(indx1,:)*D21'); % B2'*GB*[D21'; 0] (delBbigChat is 0)
    delCbigChat = real(D12'*GC(:,indx2));     % (delCbigAhat, delCbigBhat are 0)
    delCbigDhat = real(D12'*GC(:,indx1)*C2');
    delDbigDhat = real(D12'*GD*D21');         % (rest are 0)
    delAhat = delAbigAhat;
    delBhat = delAbigBhat + delBbigBhat;
    delChat = delAbigChat + delCbigChat;
    delDhat = delAbigDhat + delBbigDhat + delCbigDhat + delDbigDhat;
    g = putABCDhat(delAhat, delBhat, delChat, delDhat,structure);
else
    %USE Much more complicated calculations
    %We need Ahat, etc
    m = size(B2,2);
    p = size(C2,1);
    [Ahat, Bhat, Chat, Dhat] = getABCDhat(nhat, m, p, x,structure);
    
    F = Dhat*D22;
    %We use F = (I-Dhat*D22) a lot!
    F = eye(size(F)) - F;
    %The derivatives for Abig
    %del(f) wrt A = [G1,G2;G3,G4]  where 
    %   G1 is GA(indx1,indx1)
    %   G2 is GA(indx1,indx2)
    %   G3 is GA(indx2,indx1)
    %   G4 is GA(indx2,indx2)
    
    G1 = GA(indx1,indx1);
    G2 = GA(indx1,indx2);
    G3 = GA(indx2,indx1);
    G4 = GA(indx2,indx2);
   
    delAbigAhat = real(G4);   %This is the same for both
    delAbigBhat = real(G3*(C2+D22*(F\Dhat)*C2)'+G4*(D22*(F\Chat))');
    delAbigChat = real((B2/F)'*G2 + (Bhat*D22/F)'*G4);
    
    %This is a pretty ugly calculation and I will break it into
    %steps:    solution is lmatrix'*GA*rmatrix'
    
    lmatrix = [(B2/F); (Bhat*D22/F)];
    rmatrix = [(C2 + (D22/F)*Dhat*C2), (D22/F*Chat)];
    delAbigDhat = real(lmatrix'*GA*rmatrix');
    
    
    %Bbig
    %del(f) wrt B = [G1;G2]  where 
    %   G1 = GB(indx1,:)
    %   G2 = GB(indx2,:)
    %
    % Bbig only depends on Bhat and Dhat, others are zero
    G1 = GB(indx1,:);
    G2 = GB(indx2,:);

    %I use this quite a few times.  NO need to recompute
    D21Prime = (D21+(D22/F*Dhat)*D21);
    
    delBbigBhat = real(G2*(D21+(D22/F)*Dhat*D21)');
    delBbigDhat = real(((B2/F)'*G1 + (Bhat*D22/F)'*G2)*D21Prime');

    %Cbig 
    %del(f) wrt C = [G1 G2]  where 
    %   G1 = GB(:,indx1)
    %   G2 = GB(:,indx2)
    %
    % Cbig only depends on Chat and Dhat, others are zero
    
    G1 = GC(:,indx1);
    G2 = GC(:,indx2);
    
    %Do once and use many times
    D12Prime = D12/F;
    
    
    delCbigChat = real(D12Prime'*G2);
    %Again delCbigDhat is pretty bad, will break up into steps
    sum1 = D12Prime'*G1*(C2 + D22/F*Dhat*C2)';
    sum2 = D12Prime'*G2*(D22/F*Chat)';
    delCbigDhat = real(sum1 + sum2);
    
    %Dbig, only depends on Dhat, others 0
    delDbigDhat = real((D12/F)'*GD*D21Prime');
    
    delAhat = delAbigAhat;
    delBhat = delAbigBhat + delBbigBhat;
    delChat = delAbigChat + delCbigChat;
    delDhat = delAbigDhat + delBbigDhat + delCbigDhat + delDbigDhat;
    g = putABCDhat(delAhat, delBhat, delChat, delDhat,structure);
end
    
    
    
    
    
    
    
    