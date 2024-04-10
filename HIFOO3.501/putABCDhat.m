function x = putABCDhat(Ahat, Bhat, Chat, Dhat,structure);
% encode x from Ahat, Bhat, Chat, Dhat 

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


%IF NO STRUCTURE
if isstruct(structure) == 0
    x1 = Ahat(:);  % if use reshape must pass dimensions
    x2 = Bhat(:); 
    x3 = Chat(:);
    x4 = Dhat(:);
else
    % structure was imposed on the controller either by the user or by the
    % H2 optimization
    if (isfield(structure,'Ahat') && isfield(structure,'Bhat') && isfield(structure,'Chat')) %treat Dhat separetly
        x1 = []; x2 = [];x3 = [];
        if ~isempty(structure.Ahat)
            x1 = structure.Ahat*Ahat(:);
        end
        if ~isempty(structure.Bhat)
            x2 = structure.Bhat*Bhat(:);
        end
        if ~isempty(structure.Chat)
            x3 = structure.Chat*Chat(:);
        end
    else
        x1 = Ahat(:);  
        x2 = Bhat(:);
        x3 = Chat(:);
    end
    
    if (isfield(structure,'Dk'))
        if ~isempty(structure.Dk)        %unique solution so there is nothing to include in the optimization set
            x4 = [];
        else
            %infinity solution for Dhat
            y = Dhat(:);
            x4 = structure.V'*(y - structure.w);      % choose just the optimisation vars for H2 constraint 
        end
    elseif isfield(structure,'Dhat')                    % choose the optimization vars if no H2 constraint
        %no H2 constraint but we have Dhat structure constraint
        x4 = structure.Dhat*Dhat(:);
    else
        %no constraints on Dhat
        x4 = Dhat(:);
    end     
end

x = [x1; x2; x3; x4];
