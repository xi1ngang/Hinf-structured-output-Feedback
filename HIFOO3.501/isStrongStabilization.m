function [P,isstrongstab] = isStrongStabilization(plants)
% Run through the cellarray plants seeing if there is the 
% string "k" or "K".  
%
%Output: 
%       P : new cell array with no controller string
%       isstrongstab : 1 if strong stabilzation is selected, 0 otherwise

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

kindex =0;
isstrongstab =0;

for i = 1:numel(plant)
    curr_plant = plant{i}
    if ischar(curr_plant) && (strcmp(curr_plant,'K') || strcmp(curr_plant,'k'))
        %Found the controller
        %K cannot appear twice
        if (kindex ~= 0)
            error('Error:  Format error with input cell-array P\n')
        end
        isstrongstab = 1;
        kindex = i;
    end
end

% Now set P to be plant without the "K"


if kindex ~= 0
    P = {plant{1:(k-1)}, plant{(k+1):end}};
else
    P = plant;
end