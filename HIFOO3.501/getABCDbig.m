function [Abig, Bbig, Cbig, Dbig] = ...
    getABCDbig(A, B1, B2, C1, C2, D11, D12, D21,D22, nhat, x,structure)
% EDITED BY MARC MILLSTONE TO HANDLE D22 NOT EQUAL TO 0
% construct state space representation for general reduced-order controller
% A, B1, B2, C1, C2, D11, D12, D21 are fixed matrices
% nhat is the order of the controller
% the vector x encodes the variables Ahat, Bhat, Chat, Dhat (Ahat has
% dimension nhat by nhat, and static output feedback is the case nhat=0)
% it was originally assumed that D22 = 0 so everything is linear
% for the derivatives, see chainABCD
% from Didier's notes, Stockholm, 2005
%

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

% first decode Ahat, Bhat, Chat, Dhat from x
m = size(B2,2);
p = size(C2,1);
[Ahat, Bhat, Chat, Dhat] = getABCDhat(nhat, m, p, x,structure);

%WE need to calculate F = (I-Dhat*D22)^(-1)
% IF D22 = 0 exactly, F = I

if (D22 == 0)
       %ignore D22
        Abig = [A + B2*Dhat*C2    B2*Chat      % (n + nhat) by (n + nhat)
                Bhat*C2           Ahat    ];
        Bbig = [B1 + B2*Dhat*D21               % (n + nhat) by m1
                Bhat*D21         ];
        Cbig = [C1 + D12*Dhat*C2   D12*Chat];  % p1 by (n + nhat)   
        Dbig = D11 + D12*Dhat*D21;             % p1 by m1
else
    F = Dhat*D22;
    %We use F = (I-Dhat*D22) a lot!
    F = eye(size(F)) - F;
    Abig = [A+(B2/F)*Dhat*C2                  (B2/F)*Chat               % (n + nhat) by (n + nhat)
            Bhat*C2+(Bhat*D22/F)*Dhat*C2      Ahat+(Bhat*D22/F)*Chat];
    Bbig = [B1+(B2/F)*Dhat*D21                                       % (n + nhat) by m1
            Bhat*D21+(Bhat*D22/F)*Dhat*D21];
    Cbig = [C1 + (D12/F)*Dhat*C2   (D12/F)*Chat];  % p1 by (n + nhat)
    Dbig = D11 + (D12/F)*Dhat*D21;             % p1 by m1
end
