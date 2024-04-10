function [plant, optimize, constraints, stabilize] = processConstraintInfo(plant, obj,order,constraintInfo,options)
%This function will convert the user input into an internal structure
%useful to hifoo.  It will also perform input checking to make sure 
%that the input is correct
% If order == 0, only the objective 'h' is non-trivial for the controller, so others are pruned

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




%now verify lengths are correct
% This is partially redundant,
if length(plant) ~= length(constraintInfo)
    error('Error: The number of elements in UPPERBNDS must match the number of included plants\n')
end


if ~isempty(find(constraintInfo == -Inf))
    error('hifoo: Error: Upperbounds may not be -Inf') 
end

iscontHinf = 0;
% Now look for controller if order == 0
if order == 0
    for i = 1:numel(plant)
        thisplant = plant{i};
        if ischar(thisplant) && (strcmp(thisplant,'K') || strcmp(thisplant,'k'))
            %this is the controller
             cont_index = i;
            if (~strcmp(obj(i),'h') && ~strcmp(obj(i),'t') )
                % The controller information is implicit, so remove from everything
                if (options.prtlevel > 1)
                    fprintf('hifoo: SOF controller is implicitly stable\n')
                end
                obj = [obj(1:(i-1)) obj((i+1):end)];
                constraintInfo = [constraintInfo(1:(i-1)) constraintInfo((i+1):end)];
                plant = [plant(1:(i-1)) plant((i+1):end)];
               
            else
                iscontHinf = 1;
            end
            break;
        end
    end
end

 



optimize = find(constraintInfo == Inf); %This may be empty if no constraints requested

if ((numel(plant) == 1) && (isempty(optimize)))
    error('hifoo: If only a single plant is specified, it must be in the optimize set\n')
end



%First make sure that all plants we wish to optimize over, share a
%common constraint 

objs = obj(optimize);
if ~isempty(optimize) && max(abs(objs - circshift(objs,[0,1])))>0
    error('hifoo:  Plants that are part of the objective must have matching objectives\n')
end


%Now we set constraints.
constraints = constraintInfo(find(constraintInfo ~= Inf));

%Required to stabilize BEFORE we optimize to ensure a finite objective function
% Also, if SOF controller, we do not need to stabilize the controller
stabilize = find(obj == 'h' | obj== 'r' | obj == 't' | obj == '+');
if ((order == 0)&& (iscontHinf == 1))
    % remove controller from stabilization phase, as it is implicity stable and H-infinity norm is finite
    stabilize = stabilize(~ismember(stabilize,cont_index));
end
    


