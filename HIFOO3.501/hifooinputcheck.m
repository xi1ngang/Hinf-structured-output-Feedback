function [order,init,obj,options,upperbounds] = hifooinputcheck(inputs)
% This function processes and returns arguments of the following classes
% - order is a nonnegative integer
% - init is either a structure with fields A,B,C,D or an SS object
% - obj is a char
% - options is a structure with no required fields
% - upperbounds describes which plants belong in the objective function and which
%   belong as constraints

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

% The input arguments can be given in arbitrary order
order = [];
init = [];
obj = [];
options = [];
upperbounds = [];

for k = 1:length(inputs)
    arg = inputs{k};
       % the additional 1 is because plant is always the first arg to hifoo
    errstr = strcat('hifoo: Input argument #',int2str(k) + 1);
    switch class(arg)
    case 'double'
        %This could be either the order or the upperbounds array
        %The upperbounds array must specify at least one plant in the
        %objective function, therefore must contain at least one "Inf"
            if isempty(arg)  %ignore empty args
        
        else 
            if isempty(order)  && numel(arg) == 1
                %This is the order
                order = arg;
            elseif (isempty(upperbounds) && numel(arg) > 1)
                % more than one plant must be specified to utilize constraints
                upperbounds = arg;
            else
                error(strcat(errstr,' is not valid'));
            end
        end
     case 'struct' % [] is not a struct, it's a double, see below
            % d field is required, although a, b, c fields could be empty
            if isfield(arg, 'd')|isfield(arg,'D')|isfield(arg,'Dhat')
                if isempty(init) % i.e., init not already assigned
                    init = arg; 
                else % 2nd struct of this kind not allowed
                    error(strcat(errstr, ' is not valid'))
                end
            else
                if isempty(options)
                    options = arg;
                else
                   error(strcat(errstr, ' is not valid'))
                end
            end
    case 'char'
            if isempty(obj) % i.e., obj not already assigned
                obj = lower(arg); % convert to lower case
            else % 2nd char not allowed
                error(strcat(errstr, ' is not valid'))
            end 
			
            validchars = find(obj == 'h' | obj == 't' | obj == 'r' | obj == 's' | obj == 'p' | obj == '+');
            if length(validchars) ~= length(obj)  % this could be 1 or the number of plants 
                error(strcat(errstr, ' is a char or string but includes an invalid objective (not all h, t, r, s, or p)'))
            end
            
            % Verify that if there exists a '+' term in the objective, then
            % there is *ONLY* a '+' in the objective
            
            if sum(ismember(obj,'+')) > 0
                if(length(obj) ~= sum(ismember(obj,'+')))
                    error('The + objective cannot be combined with other objectives');
                end
            end
    case 'ss'
            if isempty(init) % i.e., init not already assigned
                init = arg; % converted to struct in hifoo
            else
                error(strcat(errstr, ' is not valid'))
            end
    otherwise
            error(strcat(errstr, ' is not valid'))
    end
    
end


%Now we set the defaults
if isempty(order)
    order = 0;
end

%We cannot set other defaults yet as that requires knowledge of the
%number of plants


