function [kinit, ssinput] = processInitialController(init,pars, options)
% Generate a controller of the proper internal format given 
% user inputted controller.  Also verifies that the given controller
% meets provided structure contraints.
% Input:  init:  The inputted controller
%         pars:  requires fields pars.m and pars.
% Outout: kinit:  The controller in the proper internal format
%         ssinput:  0 if init is structure, 1 is init is ss object.

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


ssinput = 0;
if strcmp(class(init),'ss')
    [A, B, C, D] = ssdata(init);
    init = []; % redefine init to be a structure
    init.a = A; 
    init.b = B;
    init.c = C;
    init.d = D;
    ssinput = 1; % used at end in determining output format, 
                 % overrides setting based on input "plant"
elseif isstruct(init) % default [] is not a struct
    m = pars.m;
    p = pars.p;
    if isfield(init,'a')
        % nothing to do
    elseif isfield(init,'A')
        init.a = init.A;
    elseif isfield(init,'Ahat')
        init.a = init.Ahat;
    else
        init.a = zeros(0,0); % an order 0 (SOF) controller
    end
    initorder = size(init.a, 1);
    if initorder ~= size(init.a, 2)
        error('hifoo: input structure "init.a" must be square')
    end
    if isfield(init,'b')
        % nothing to do
    elseif isfield(init,'B')
        init.b = init.B;
    elseif isfield(init, 'Bhat')
        init.b = init.Bhat;
    else
        init.b = zeros(0,p); % an order 0 (SOF) controller
    end
    if size(init.b, 2) ~= p
        error('hifoo: input "init.b" must have column dimension %d', p)
    elseif size(init.b, 1) ~= initorder
        error('hifoo: inputs "init.a" and "init.b" must have same row dimension')
    end
    if isfield(init,'c')
        % nothing to do
    elseif isfield(init,'C')
        init.c = init.C;
    elseif isfield(init, 'Chat')
        init.c = init.Chat;
    else
        init.c = zeros(m,0); % an order 0 (SOF) controller
    end
    if size(init.c, 1) ~= m
        error('hifoo: input "init.c" must have row dimension %d', m)
    elseif size(init.c, 2) ~= initorder
        error('hifoo: inputs "init.a" and "init.c" must have same column dimension')
    end
    if isfield(init,'d')
        % nothing to do
    elseif isfield(init,'D')
        init.d = init.D;
    elseif isfield(init,'Dhat')
        init.d = init.Dhat;
    else % this is the only required field as it cannot be empty
        error('hifoo: input structure "init" must have a field "d"')
    end
    if any(size(init.d) ~= [m p]) % for example, if init.d is empty
        error('hifoo: input "init.d" must have dimension %d by %d', m, p)
    end
    ssinput = 0; % used at end in determining output format, 
                 % overrides setting based on input "plant"
end % if init == [], ssinput value is set depending on input "plant"


% Verify that the inputted controller as the proper structure
if isstruct(init) && isstruct(options.structure)
    if (isfield(options.structure,'Ahat') && isfield(options.structure,'Bhat')...
            && isfield(options.structure,'Dhat') && isfield(options.structure,'Chat') )
        %doing weird transformations because of the new way the structure
        %is treated
        if ~isempty(options.structure.Ahat)
            A = options.structure.Ahat'*ones(size(options.structure.Ahat,1),1);
        else
            A = ones(size(init.a));
        end
        
        if ~isempty(options.structure.Bhat)
            B = options.structure.Bhat'*ones(size(options.structure.Bhat,1),1);
        else
            B = ones(size(init.b));
        end
        
        if ~isempty(options.structure.Chat)
            C = options.structure.Chat'*ones(size(options.structure.Chat,1),1);
        else
            C = ones(size(init.c));
        end
        
        if ~isempty(options.structure.Dhat)
            D = options.structure.Dhat'*ones(size(options.structure.Dhat,1),1);
        else
            D = ones(size(init.d));        %dimensions were already checked
        end
        
        % check if the initial controller is the same order as the
        % structure
        
        if (sum(size(A) ~= size(init.a)) > 0)
            error('The given initial controller doesn''t have the same order as the given structure');
        end
        
        % check element by element if the fixed parameters are respected in
        % the initial controller - works for old and new structure
	
        r = not(logical(reshape(A,size(init.a))));
        if max(abs(init.a(r) - options.structure.Afix(r))) > 0
            error('The initial controller matrix a differs from the user-provided structure');
        end
        
        r = not(logical(reshape(B,size(init.b))));
        if max(abs(init.b(r) - options.structure.Bfix(r))) > 0
            error('The initial controller matrix b differs from the user-provided structure');
        end
        
        r = not(logical(reshape(C,size(init.c))));
        if max(abs(init.c(r) - options.structure.Cfix(r))) > 0
            error('The initial controller matrix c differs from the user-provided structure');
        end
        
        r = not(logical(reshape(D,size(init.d))));
        if max(abs(init.d(r) - options.structure.Dfix(r))) > 0
            error('The initial controller matrix d differs from the user-provided structure');
        end

	%Already known to be the proper format
%        if nnz(logical(init.a) - logical(reshape(A,size(init.a)))) > 0 || ...
%                nnz(logical(init.b) - logical(reshape(B,size(init.b)))) > 0 || ...
%                nnz(logical(init.c) - logical(reshape(C,size(init.c)))) > 0 || ...
%                nnz(logical(init.d) - logical(reshape(D,size(init.d)))) > 0
%            
%            error('The user-provided structure differs from the given initial controller');
%        end
    end
end

kinit = init;
