function verifyStabilityBoundary(K,pars,prtlevel)
% Verify that the eigenvalues of the Abig
% matrix are remaining within the stability boundary
% Let x+yi be the eigenvalues of Abig.  We will
% ensure that x/(1+x^2+y^2) <= 1e-6

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

structure = pars.structure;
nhat = pars.nhat;
warn = 0;
for i = 1:numel(pars.plantinfo)
    plant = pars.plantinfo{i};
    if (isfield(plant,'controller') && plant.controller == 1)
           %This is for strong stabilization, this is the controller

           m = pars.m;
           p = pars.p;

           [Abig, Bbig, Cbig, Dbig] = getABCDhat(pars.nhat, m, p, K,pars.structure);
    else

       A = plant.A; B1 = plant.B1; B2 = plant.B2; C1 = plant.C1; C2 =plant.C2; 
       D11 = plant.D11; D12 = plant.D12; D21 = plant.D21; D22 = plant.D22;
       structure = pars.structure;
       n = size(A,1);

       [Abig, Bbig, Cbig, Dbig] = ...
       getABCDbig(A, B1, B2, C1, C2, D11, D12, D21,D22, nhat, K,structure );
    end
       
    evals = eig(Abig);
    boundary = real(evals)./(1+real(evals).^2+imag(evals).^2);
    if nnz(boundary >= -1e-6) > 0
        warn = warn + 1;
    end
end
if warn > 0 && prtlevel > 0
    fprintf('hifoo: An eigenvalue of the closed-loop system is approaching the stability boundary\n');
end