function [x,N,x0] = solveSystem(A,b)
%   solves the linear sistem Ax = b
%   returns x = nan if no solution, x = [] for infinity solution or x if
%   unique solution
%   for the infinity solution x = x0 + N*y where y is arbitrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HIFOO, A Matlab package for Fixed Order H-infinity and H2 control 
%%  Copyright (C) 2010  Georgia Deaconu
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

[m,n] = size(A);
[m,p] = size(b);

[U,S,V] = svd(A);   %   A = U*S*V'

%solve the system S*y = c with y = V'*x and c = U'*b

c = U'*b;
%r = rank(S);
r = rank(S,1e-10*S(1,1));

idx = [];
if r < m                    %then we should check the norm in c
    idx = r+1:m;
end

if norm(c(idx)) > 1e-10*S(1,1)        %no solution - threshold chosen randomly
    x = nan;
    N = nan;
    x0 = nan;
    nf = nan;
    return;
else
    nf = n-r;
    y = zeros(r,p);
    % unique solution or infinity solution
    for i = 1:r
        y(i) = c(i)/S(i,i);
    end
        
    if (nf == 0)            %unique solution
        x = V*y;
        x0 = [];
        N = [];
    else
        %infinity solution
        x0 = V*[y; zeros(nf,1)];
        x = [];
        N = V*[zeros(r,nf); eye(nf)];
    end
end
