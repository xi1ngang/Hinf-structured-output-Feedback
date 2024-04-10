function xinit = xinitsetup(m, p, order, init, obj, pars,structure)
% Intended to be called by hifoomain only.
% Set up initial point for optimization. 
% init is the initial controller, "order" is the desired controller order.
% First code init as a vector xinit, with order initorder. Then, if needed,
% expand or truncate xinit.  If order > initorder (the "expand" case)
% the expansion of xinit corresponds to padding Bhat and Chat with zeros
% and padding Ahat with zeros too except that its additional diagonal 
% entries are set to a negative number, say diagpad. This ensures
% that Abig is also padded with zeros except for the same additional
% diagonal entries, diagpad. In fact, we set diagpad = diagpad0 - shift,
% where shift is an arbitrarily chosen positive number. If the
% current Abig is itself unstable, we can just set diagpad0 = 0.
% If the current Abig is stable, the choice of diagpad0 depends on the
% objective.  We can set it to the current objective value for Abig
% in the case of the spectral and pseudospectral abscissa.  In the
% case of the H-infinity norm objective, anything is OK, since it will
% be annihilated by the padded zeros in Bhat and Chat.  However, in
% the case of the inverse-of-complex-stability-radius objective, we
% need to set diagpad0 to minus the complex stability radius of the
% current Abig (not the inverse; the reason we chose the stability 
% radius objective to be the inverse, as opposed to minus, the radius, 
% was simply to be consistent with the H-infinity norm defintion). 
% This ensures that none of the optimization objectives is increased 
% by the increase in order. 
%
% Note: if there is a nonzero penalty term it DOES increase 
% since A is padded with a nonzero diagonal. If there is a barrier
% term (a multiple of the inverse stability radius added to an 
% H-infinity objective) we can ensure that this combined objective does 
% not increase by treating this case as if the objective is the inverse
% stability radius, since, as noted above, the new diagonal entries don't 
% affect the H-infinity objective as they are annihilated by the padded 
% zeros in Bhat and Chat.
% 
% If order < initorder, x is truncated arbitrarily; in this case 
% the optimization objectives change arbitrarily.  
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



% hifoo already checked that dimensions are compatible
xinit = putABCDhat(init.a, init.b, init.c, init.d,structure);
initorder = size(init.a,1); % the case 0 is an SOF initial controller

if initorder ~= order  % if same, nothing else to do
    % if we have a structure imposed on the controller and this structure
    % is not imposed by the H2 optimization
    if isstruct(structure) && ~isfield(structure,'Dk') 
            error('Error: When using structure, initial controller order must match desired controller order');
    else
    
    incr = order - initorder;  
    % Expand (or contract) xinit.
    % If initorder < order, none of the optimization objectives change.
    % If initorder > order, optimization objective may change arbitrarily.
    % Set Ahat, ...,Dhat to init.a,...init.d
    [Ahat, Bhat, Chat, Dhat] = getABCDhat(initorder, m, p, xinit,structure);
    if incr > 0
        Bhat = [Bhat; zeros(incr,p)];  % append zero row(s) to Bhat
        Chat = [Chat zeros(m,incr)];   % append zero col(s) to Chat
        % pad Ahat with diagonal entry (entries)
        shift = 1; % arbitrary positive number
        % alternatives, which are not defined if initorder is 0, are 
        % shift = min(1, max(max(abs(Ahat)))); % may be bad if Ahat is 0 
        % shift = max(1, max(max(abs(Ahat)))); % may be bad if Ahat is big 
        %
        % the following changes to pars are not passed back to hifoomain
        % compute spectral abscissa of current Abig (defined by current xinit)
        pars.objective = 's';
        pars.nhat = initorder;
        pars.penalty = 0; % see note above; pars.barrier will not be accessed
        f0 = hifooobj(xinit, pars); 
        if f0 >= 0 % happens if and only if current Abig is unstable
            diagpad0 = 0; % avoid introducing a new positive eigenvalue for Abig
        elseif obj == 'r' | obj == 'h'  
            % compute inverse of complex stability radius of current Abig
            pars.objective = 'r';
            f0 = hifooobj(xinit, pars); 
            diagpad0 = -1/f0; % negative of stability radius, small if nearly unstable
        else
            diagpad0 = f0; % must be negative; see long comment at start of code
        end
        for j=1:incr
            Ahat(initorder + j, initorder + j) = diagpad0 - shift;
        end
    else
        indx = 1:order;
        Ahat = Ahat(indx,indx);  % arbitrary truncations
        Bhat = Bhat(indx,:); 
        Chat = Chat(:,indx);
    end
    end
    xinit = putABCDhat(Ahat, Bhat, Chat, Dhat,structure);
end