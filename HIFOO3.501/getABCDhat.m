function [Ahat, Bhat, Chat, Dhat, varsA, varsB,varsC,varsD] = getABCDhat(nhat, m, p, x,structure);
% decode Ahat, Bhat, Chat, Dhat from x

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
if isstruct(structure) == 0
    %fprintf('getABCDhat: we don''t have any structure \n');
    indx1 = 1:nhat^2;
    indx2 = nhat^2 + 1 : nhat^2 + nhat*p;
    indx3 = nhat^2 + nhat*p + 1 : nhat^2 + nhat*p + nhat*m;
    indx4 = nhat^2 + nhat*p + nhat*m + 1 : nhat^2 + nhat*p + nhat*m + m*p;
    if length(x) ~= indx4(length(indx4))
        fprintf('getABCDhat: length mismatch')
    end
    Ahat = reshape(x(indx1), nhat, nhat);
    Bhat = reshape(x(indx2), nhat, p);
    Chat = reshape(x(indx3), m, nhat);
    Dhat = reshape(x(indx4), m, p);
    
    varsA = nhat*nhat;
    varsB = nhat*p;
    varsC = m*nhat;
    varsD = m*p;
else
    %WE ADD STRUCTURE GIVEN BY structure struct
    if isfield(structure,'Ahat') && isfield(structure,'Bhat') && isfield(structure,'Chat')
        %Structure for Ahat, Bhat, Chat 
        varsA = size(structure.Ahat,1);
        varsB = size(structure.Bhat,1);
        varsC = size(structure.Chat,1);
        
        indx1 = 1:varsA;
        indx2 = varsA+1 : varsA+varsB;
        indx3 = varsA+varsB+1 : varsA+varsB+varsC;
        
        %Multiply by the transpose to put the zeros in the correct position for the 
        % elements that do not enter the optimization
        if ~isempty(structure.Ahat)
            Ahat = reshape(structure.Ahat'*x(indx1),nhat,nhat) + structure.Afix;
        else
            Ahat = zeros(nhat,nhat) + structure.Afix;
        end
        if ~isempty(structure.Bhat)
            Bhat = reshape(structure.Bhat'*x(indx2),nhat,p) + structure.Bfix;
        else
            Bhat = zeros(nhat,p) + structure.Bfix;
        end
        if ~isempty(structure.Chat)
            Chat = reshape(structure.Chat'*x(indx3),m,nhat) + structure.Cfix;
        else
            Chat = zeros(m,nhat) + structure.Cfix;
        end
    else
        %no structure so we just reshape everything
        indx1 = 1:nhat^2;
        indx2 = nhat^2 + 1 : nhat^2 + nhat*p;
        indx3 = nhat^2 + nhat*p + 1 : nhat^2 + nhat*p + nhat*m;
        
        Ahat = reshape(x(indx1), nhat, nhat);
        Bhat = reshape(x(indx2), nhat, p);
        Chat = reshape(x(indx3), m, nhat);
        
        varsA = nhat*nhat;
        varsB = nhat*p;
        varsC = m*nhat;
    end
    
    % if Dk exists it includes the information from Dhat
    if isfield(structure,'Dk')
        if ~isempty(structure.Dk)
            % unique solution for Dhat
            Dhat = structure.Dk;
        else
            % infinity solution for Dhat
            indx4 = varsA + varsB + varsC + 1:varsA + varsB + varsC + size(structure.V,2);
            y = x(indx4);
            d = structure.V*y + structure.w;
            Dhat = reshape(d,m,p);
        end
    elseif isfield(structure,'Dhat')
        % no H2 constraint, structure imposed by user for Dhat 
        varsD = size(structure.Dhat,1);
        indx4 = varsA+varsB+varsC+1: varsA+varsB+varsC+varsD;
        
        if(length(x) ~= varsA+varsB+varsC+varsD )
            fprintf('getABCDhat:length mismatch with structure\n');
        end
        
        Dhat = reshape(structure.Dhat'*x(indx4),m,p) + structure.Dfix;
    else
        %no constraints so just reshape
        indx4 = nhat^2 + nhat*p + nhat*m + 1 : nhat^2 + nhat*p + nhat*m + m*p;
        if length(x) ~= indx4(length(indx4))
            fprintf('getABCDhat: length mismatch')
        end
        Dhat = reshape(x(indx4), m, p);
        varsD = m*p;
    end
    
end
