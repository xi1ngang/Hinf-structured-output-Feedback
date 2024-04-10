function [K, f, viol, loc] = hifoo(P, varargin)
%
%  HIFOO, A Matlab package for Fixed Order H-infinity and H2 Control 
%   Stabilization and Performance Optimization for Multiple Plants
%
%   K = HIFOO(P) looks for a static output feedback controller that
%    stabilizes the plants P{j},j=1,2,... and locally minimizes the max of
%    the H-infinity norms of the closed-loop plants.
%   K = HIFOO(P, ORDER) looks for an order ORDER controller doing the same.
%   K = HIFOO(P, 't') looks for a static output feedback controller that
%    stabilizes the plants P{j},j=1,2,... and locally minimizes the max of
%    the H2 norms of the closed-loop plants.
%   K = HIFOO(P, ORDER, 't') or K = HIFOO(P, 't', ORDER) looks for an
%    order ORDER controller doing the same.
%   K = HIFOO(P, 'r') looks for a static output feedback controller that
%    stabilizes the plants P{j},j=1,2,... and locally minimizes the max of
%    reciprocals of the complex stability radii of the closed-loop plants.
%   K = HIFOO(P, ORDER, 'r') or K = HIFOO(P, 'r', ORDER) looks for an 
%    order ORDER controller doing the same.
%   K = HIFOO(P, 's') looks for a static output feedback controller that 
%    locally minimizes the max of the spectral abscissas of the closed-loop 
%    plants (max of the real parts of their eigenvalues).
%   K = HIFOO(P, ORDER, 's') or K = HIFOO(P, 's', ORDER) looks for an
%    order ORDER controller doing the same.
%   [K, F, VIOL, LOC] = HIFOO(P, ORDER, INIT, FUN, UPPERBND, OPTIONS)
%    looks for an order ORDER controller for P{j},j=1,2,... that 
%    locally solves the following optimization problem:
%		min F(K) = max{g_j(K):  UPPERBND(j) = +inf}
%        subject to the constraints 
%             g_j(K) <= UPPERBND(j),  j=1,2,...,
%    where each g_j is one of the supported closed-loop functions 
%    encoded in FUN, using an initial guess INIT, with options 
%    given in OPTIONS, returning the controller K, objective value F, 
%    local optimality certificate LOC and constraint violations VIOL. 
%
%  Input Parameters:
%   P is a cell array of plants, each of which has one of 4 formats
%     (if there is only one plant, either P or P{1} may be set to one
%      of the first 3 formats):
%     - structure with fields 
%        A,B1,B2,C1,C2,D11,D12,D21,D22
%         B can be used instead of B2, C can be used instead of C2,
%         D can be used instead of D22.
%         If D or D22 is not provided, it is assumed to be zero. 
%        (B1,C1,D11,D12,D21 are required only if the supported function
%         corresponding to the plant is H-infinity performance, and 
%         are ignored otherwise, OR
%     - SS object (as defined in MATLAB's Control System Toolbox), 
%        encoding A,[B1 B2],[C1; C2],[D11 D12; D21 D22] in standard format 
%        (for computing H-infinity performance, the index partitioning 
%         of [B1 B2] and [C1 C2] should be specified in InputGroup.U1, 
%         InputGroup.U2, OutputGroup.Y1 and OutputGroup.Y2, which are
%         set and viewed by "set" and "get" respectively), OR
%     - string giving the COMPleib name of the plant 
%        (see www.compleib.de), OR
%     - the single character 'K', equivalent to providing a plant for
%         which the closed loop plant equals the controller K.
%   ORDER, INIT, FUN, UPPERBND and OPTIONS are optional and may be
%     specified in any permutation order
%    ORDER: the order of the controller, a nonnegative integer
%     (default: 0 (static output feedback))
%    INIT: initial guess for controller, has one of 2 formats:
%     - structure with fields: a, b, c, d, OR   
%     - SS object
%     When OPTIONS.struct is not used, if the order of the initial guess is
%     less than the desired order, the initial guess is augmented to have the 
%     desired order without increasing the objective value. Thus this routine 
%     can be called repeatedly to get successively better controllers as the
%     order is increased. If the order of the initial guess is greater than 
%     the desired order, the initial guess is truncated to have the desired
%     order but this may increase the objective value.
%     If OPTIONS.struct is used, the order of the initial controller must 
%     equal ORDER.
%     (default: INIT is generated randomly)
%    FUN: a character or a string, specifying the functions g_j.
%      If FUN is a single character, it specifies one supported function
%      from the following list that applies to all plants.
%      If FUN is a string, it must have length equal to the number of
%      plants and it must specify a supported function from the following
%      list for each plant; however, plants with UPPERBND(j) = +inf must
%      all be associated with the same supported function.
%      The supported functions are:
%       'h': +inf if closed loop plant is unstable, otherwise
%            H-infinity norm of transfer function from performance input 
%            to performance output for the closed loop plant
%       't': +inf if closed loop plant is unstable, otherwise H2 norm of
%            transfer function from performance input to performance output
%            for the closed loop plant
%       'r': +inf if closed loop plant is unstable, otherwise
%            reciprocal of complex stability radius of closed loop plant
%       's': spectral abscissa (max(real(eigenvalues))) of closed loop plant
%       'p': pseudospectral abscissa of closed loop plant 
%            (a value for epsilon must be provided in OPTIONS.epsilon)     
%       '+': similar to 's', except that spectral abscissa is not optimized:
%            hifoo only attempts to find a controller that simultaneously 
%            stabilizes *all* closed loop plants. Consider using 's' instead.
%      (default: 'h')
%    UPPERBND: a vector specifying upper bounds on the g_j. 
%       All plants for which UPPERBND(j) = +inf are unconstrained and
%       enter the minimization objective instead.
%       Note that if FUN(j) is 'h', 't' or 'r', a stability constraint on
%       the closed loop plant for P{j} is implicit as the value of g_j is
%       +inf if the closed loop plant is unstable. 
%       Note: it is allowable to constrain all plants, that is set all
%       components of UPPERBND to be finite, unless there is only one plant. 
%      (default: all +inf)
%    OPTIONS: structure with all fields optional:
%       OPTIONS.cpumax:  quit when cpu time in seconds exceeds this
%          (must be a positive number; default: inf)
%       OPTIONS.fast: 
%          1 to use a fast optimization method (BFGS) only 
%          0 to finish optimization with slower method (gradient sampling,
%            which may give a better answer and a better output LOC)
%          (default: 0)
%       OPTIONS.prtlevel: one of 0 (no printing), 1 (minimal, default),
%          2 (includes output from HANSO), 3 (verbose)
%       OPTIONS.struct:  a structure of matrices with fields
%          OPTIONS.struct.a, OPTIONS.struct.b, 
%          OPTIONS.struct.c, OPTIONS.struct.d that can be used
%          constrain individual entries of the controller matrices
%          K.a, K.b, K.c, K.d to have specified values.
%          In previous versions of HIFOO, these entries must all be either
%          0 or 1, with 0 indicating that the corresponding value in the
%          controller matrix is fixed at 0 and 1 indicating that the value 
%          is free to vary.  In version 3.5 and higher, any fixed value 
%          can be specified by setting the corresponding entry to that
%          value, and variables that are free are indicated by setting the 
%          corresponding entries to NaN.  Both conventions are supported
%          by versions 3.5 and higher.  Which convention is being used
%          is recognized by checking for NaNs in OPTIONS.struct.a,...
%          OPTIONS.struct.d. If there are none, then 1's indicate free
%          values and 0's indicate fixed values.  If there are NaN's, then
%          these indicate free values and all other numbers indicate
%          fixed values.
%       OPTIONS.nrand: number of starting points to use in addition to INIT
%          (random perturbations of increasing sizes if INIT is provided)
%          (must be a positive integer; default: 3)
%       OPTIONS.normtol: termination tolerance (see output parameter LOC)
%          (default: 1e-3)
%       OPTIONS.evaldist: evaluation distance (see output parameter LOC)
%          (default: 1e-3)
%       OPTIONS.weight: a vector of positive weights, one for each plant.
%          If w = OPTIONS.weight is provided, the problem solved is
%		      min F(K) = max{w(j)g_j(K):  UPPERBND(j) = +inf}
%              subject to the constraints 
%                   w(j)*g_j(K) <= UPPERBND(j),  j=1,2,...,
%          (default: all ones) 
%       OPTIONS.weightNormK:  weight for adding a penalty on the size of
%          the controller to the objective function, specifically 
%          sqrt(||K.a||^2 + ||K.b||^2 +  ||K.c||^2 + ||K.d||^2),  
%          where K.a, ..., K.d are the controller matrices
%          (default: 0)
%       OPTIONS.augmentHinf: weight for adding reciprocal of complex
%          stability radius to H-infinity norm function to avoid closed 
%          loop plants that are only marginally stabilized: applies to all 
%          plants for which g_j is the H-infinity norm. 
%          (default: 0)
%       OPTIONS.epsilon: in case that any FUN(j) is 'p', value of epsilon
%          defining epsilon-pseudospectral abscissa
%          (default: 0.01)
%    
%  Output parameters
%   K: best controller found, has one of 2 formats:
%     - structure with 4 fields: a, b, c, d
%     - SS object encoding a, b, c, d in SS format
%     the format is compatible with the format used for INIT
%     (or the format used by P if INIT is not provided)
%   F: corresponding value of minimization objective (including any 
%      terms arising from options.weightNormK and options.augmentHinf,
%      but not including any penalty for constraint violations)
%   VIOL: a vector specifying constraint violations for the 
%      plants: VIOL(j) = max(0, g_j(K) - UPPERBND(j))
%      (or max(0, w(j)*g_j(K) - UPPERBND(j)) if not unit weights w(j) are
%      provided in OPTIONS.weight)
%   LOC: local optimality certificate, structure with 2 fields, both set
%      to NaN if F is INF or VIOL is nonzero.  Otherwise:
%       LOC.dnorm: norm of a vector in the convex hull of bundled or
%          sampled gradients of the minimization objective,
%          evaluated at and near K 
%       LOC.evaldist: specifies max distance from K at which these 
%          gradients were evaluated.  If this is 0, then LOC.dnorm is 
%          just the norm of the gradient of the objective at K.
%      The smaller LOC.dnorm and LOC.evaldist are, the more likely 
%       it is that K is an approximate local minimizer.
%
%  *Useful Tip*: if the output K is not satisfactory, try calling HIFOO 
%    again with INIT set to K.
%
%   Other software needed
%    Required: HANSO 2.0 (Hybrid Algorithm for Non-Smooth Optimization),
%      www.cs.nyu.edu/overton/software/hanso/hanso2_0/
%    Required when any g_j is 'h' or 'r': either MATLAB's Control System 
%      Toolbox or SLICOT's linorm (see options.hinfalg, below)
%    Required when any g_j is 't': MATLAB's Control System Toolbox  
%
% Version 3.5, 2011, see GPL license info below.
% Further documentation is at www.cs.nyu.edu/overton/software/hifoo/
% Marc Millstone, Michael L. Overton, Didier Henrion, Georgia Deaconu,
%  Suat Gumussoy, Denis Arzelier
% Send comments/bug reports to Marc Millstone and Michael Overton,
% lastname@cims.nyu.edu, with a subject header containing "hifoo". 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HIFOO, A Matlab package for Fixed Order H-infinity and H2 control 
%%  Copyright (C) 2010  Marc Millstone, Michael Overton, Georgia Deaconu
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% end of help file%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  Additional options which are not intended for typical users, and 
%  therefore not part of the help info.  They may be useful to us or to
%  advanced users.
%       OPTIONS.rho: a structure used to define the parameter rho
%          in the L1 penalty function used to impose the constraints,
%          with 3 fields (see below for more information)
%          OPTIONS.rho.init: initial value for rho 
%          OPTIONS.rho.multiply: the amount rho is 
%           multiplied each time a larger value is tried
%          OPTIONS.rho.max: The maximum value allowed for
%           rho before quitting
%          (default: 100, 100 and 1e10 respectively)
%       OPTIONS.maxit:  The maximum number of iterations allowed for BFGS
%       OPTIONS.hinftol:  the tolerance used by the H-infinity norm 
%           computation (default: 1e-7)
%       OPTIONS.hinfalg:  a string which selects the algorithm for
%            computing the H-infinity norm
%           'matlab' -- Uses the H-infinity norm provided in the Matlab 
%                       Robust Control Toolbox.  Experiments show that this
%                       is slower, but more robust
%           'slicot' -- Uses H-infinity norm, linorm, provided by SLICOT 
%                       (free for non-commercial use).  To use this 
%                       functionality, linorm must be in the path.
%                       Experiments show that this is faster, however,
%                       a small error may be observed.
%           (default ='matlab')
%%%%%%%%%%%%%%%%%%%%%%%%%%%Information about the Method %%%%%%%%%%%%%%%%%%
%     There are two phases, Stabilization and Optimization
%   STABILIZATION PHASE
%     The purpose of this phase is to stablize all plants P{j} for which
%   FUN(j) is 'h', 'r' or 't' (H-infinity norm or reciprocal of complex 
%   stability radius or H2 norm, respectively).  Stabilizing these plants, 
%   we obtain an initial point for which the corresponding g_j is FINITE. 
%   We apply a combination of BFGS and Gradient Sampling to the max of the
%   spectral abscissas (max(real(eig))) of the corresponding closed loop 
%   plants, starting from INIT and OPTIONS.nrand successively larger random 
%   perturbations of INIT, or OPTIONS.nrand random starts if INIT is not 
%   provided, until finding one or more controllers that stabilize these
%   plants, or default iteration limits are exceeded.   
%   If this phase fails, HIFOO terminates with an error message and with 
%   F or some VIOL(j) (or both) equal to +inf. The stabilization 
%   phase is skipped if no plants are associated with 'h', 'r' or 't'.
%   OPTIMIZATION PHASE
%      Using the starting points found by the stabilization phase, or,
%   if this was skipped, INIT and OPTIONS.nrand random perturbations,
%   we try to optimize the L1 penalty function
%         max{g_j(K): UPPERBND(j) = inf} + rho*sum{VIOL(j)}
%   where VIOL(j) = max(0, g_j(K) - UPPERBND(j)).
%   (If options.weight is not all 1s, weights are included as above;
%   if options.weightNormK is not 0, there is an additional term.)
%   The initial value for rho is OPTIONS.rho.init.  The value of the L1
%   penalty function is finite at the starting points and points generated 
%   with infinite objective value will be rejected by the line search. 
%   After the minimization is done, a check is made as to whether any 
%   VIOL(j) is nonzero, and if so, rho is multiplied by OPTIONS.rho.multiply, 
%   and we try again.  This is repeated if necessary until a feasible point
%   is found or rho is larger than the maximum value allowed, 
%   OPTIONS.rho.max, when termination takes place, returning  VIOL(j).  
%   We use only BFGS (not gradient sampling) until we find a feasible point, 
%   and then call HANSO with the final value for rho to try to 
%   improve this using gradient sampling (this might possibly require
%   increasing the penalty parameter some more.)
%   (There is no point in repeating with a larger rho if there are
%   no terms in the objective function, that is ALL plants are constrained.)
%
%   HANSO 2.0, Hybrid Algorithm for Non-Smooth Optimization, is a
%   hybrid of BFGS and gradient sampling. Termination takes place when
%   LOC.dnorm <= OPTIONS.normtol with LOC.evaldist approximately <= 
%   OPTIONS.evaldist, or default iteration limits are exceeded, or 
%   OPTIONS.cpumax CPU time is exceeded. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF MAIN COMMENTS%%%%%%%%%%%%%%%%%%%%%%%%

% first decode the optional arguments provided to varargin
[order, init, obj, options,constraintInfo] = hifooinputcheck(varargin);
if ~isfield(options, 'prtlevel')
    options.prtlevel = 1; % other fields checked in hifoomain
end 
prtlevel = options.prtlevel; 
if prtlevel > 0 % options fields printed by hifoomain
    fprintf('HIFOO Version 3.5, using options\n')
end
% convert input format for plant, making all the necessary input checks
ssinput = 0;
%See if the inputted plant is a plant or a cell array
% If the former, convert to cell array.

if (~iscell(P))
    P = {P};
end
pars.nplants = numel(P);

%If no constraint information, then everything is in 
% the optimize variable
if isempty(constraintInfo)
    constraintInfo = Inf*ones(1,numel(P));
end

if isempty(obj)  
    %The default objective is 'h'
    obj = char('h'*ones(1,numel(P)));
end

if (length(obj)== 1  &&  pars.nplants > 1) 
    % User utilized shortcut for calling hifoo
    obj = char(obj*ones(1,pars.nplants));
elseif length(obj) ~= pars.nplants
    error('hifoo: length of obj string is not 1 or number of plants')
end
    

%Set indices for which plants belong in objective and which are constraints
% verify proper formatting
[P, options.optimize,options.constraints, options.stabilize] = processConstraintInfo(P,obj,...
        order, constraintInfo,options);


% Process the plant input.  Also verify that the dimensions are compatible
[pars.plantinfo, isstrongstab, ssinputplant, pars.m,pars.p] = processPlants(P,constraintInfo,obj,options);



% Set the structure for the controller
[options.structure] = processStructure(pars,options,order);

% now process the init argument 
% if it was omitted, it was already set to [] by hifooinputcheck, and
% this will be checked by hifoomain
% Verifies that the initial controller has the proper structure
[init,ssinputcont] = processInitialController(init,pars,options);


%If user desires strong stabilization, add controller as a plant
%THIS TECHNIOQUE IS DEPRECATED USE "K" as a plant instead
% Also pars.strongstabweight is deprecated, use options.weights() instead
if (isfield(options,'strongstab'))
    
    if options.prtlevel > 0
        warning('UPDATE: This technique for using strong stabilization is deprecated. See comments\n') 
    end
    if options.strongstab == 1
        %Need to stabilize controller as well.  Taken care of in 
        %hifooobj
        options.optimize = [options.optimize, options.optimize(end)+1];
        pars.plantinfo{numel(pars.plantinfo)+1} = -1;
        if (isfield(options,'strongstabweight'))
            pars.strongstabweight = options.strongstabweight;
        else
            pars.strongstabweight = 1;
        end
    end
end



if isstrongstab == 1
      %Need to stabilize controller as well.  set up appropriately above
        if (isfield(options,'strongStabWeight'))
            pars.strongstabweight = options.strongstabweight;  % CHANGE THIS TO MATCH
        else
            pars.strongstabweight = 1;
        end
end
    

%Set objective to be whatever the optimize objective is
% first make sure it is not empty
if ~isempty(options.optimize)
    obj = obj(options.optimize(1));
else
    %Everything is a constraint, so set it to be whatever
    obj = 'c';
end


% If there is only one plant and no strong stabilization,
% the maximum order required is that of the plant itself.

if (numel(pars.plantinfo)  == 1 && ~isstrongstab)
    n = length(pars.plantinfo{1}.A);
    if order > n
      error(sprintf('hifoo: order %d is greater than state space dimension %d\n', order, n))
    end
end

% identify the case of H2 SOF controller with unique solution
if  isfield(options.structure,'V')
    %     we have h2 objective
    if isfield(options.structure,'D22') &&  (sum(sum(options.structure.D22 ~= zeros(size(options.structure.D22)))) ~= 0)
        % we are in the case where we have h2 objective and D22
        % it is not possible to have also a defined structure for the controler 
            if isfield(options.structure,'Dhat')
                error('hifoo: cannot have H2 objective, D22 != 0 and a structure constraint on the controller\n');
            end
        end
        
        if isnan(options.structure.Dk)
            error('hifoo: the condition for the existence of the finite H2 controler is not satisfied\n');
        elseif ~isempty(options.structure.Dk)  % unique solution for Dk
            
            M = [];
            
            if (sum(sum(options.structure.D22 ~= zeros(size(options.structure.D22)))) ~= 0)  % D22 != 0
                N = options.structure.Dk*options.structure.D22;
                I = eye(size(N));
                M = I + N;     % 1+Dk*D22
                
                cd = cond(M);      % the conditionment number for the patrix that should be inverted
                
                if cd*eps>1
                    if cd == Inf    % 1+Dk*D22 == 0
                        error('hifoo: unique solution for Dk and the transformation matrix 1+Dk*D22 == 0, one of the closed-loop plants may not be well-posed ');
                    else
                        warning('hifoo: unique solution for Dk and the transformation matrix 1+Dk*D22 is singular to working precision, one of the closed-loop plants may not be well-posed \n');
                    end
                end
            end
            
            if (order == 0)       %H2: unique solution for Dk and SOF controller, nothing else to do
                fprintf('hifoo: unique solution for the SOF controller \n');
                
                pl = pars.plantinfo{1};
                m = size(pl.B2,2);
                p = size(pl.C2,1);
                
                if isempty(M)
                    K.d = options.structure.Dk;
                else
                    K.d = M\options.structure.Dk;
                end
                
                K.a = zeros(order,order); K.b = zeros(order,p); K.c = zeros(m,order);
                
                if (ssinputplant) % convert to ss object format of Control Systems Toolbox
                    K = ss(K.a, K.b, K.c, K.d);
                end
                
                % compute the objective function using hifooobj
                x = K.d(:);
                
                pars.slicot = 0;    % we don't need it but hifooobj wants it
                pars.nhat = order;
                [pars,options] = processOptions(pars,options);
                pars.penalty  = options.penalty;
                              
                pars.objective = 'i'; % i here means independent or individual.  hifooobj will look at the set objective in pars.plantinfo
                %  No one objective function, each plant carries its own
                
                % We compute the objective value with only terms in the
                % actual objective portion of the problem.
                parsToEvaluateObjective = pars;
                parsToEvaluateObjective.nplants = length(pars.optimize);
                parsToEvaluateObjective.plantinfo = {pars.plantinfo{pars.optimize}};
                
                f = hifooobj(x,parsToEvaluateObjective)
                return
            end
        end
    end



%solve the problem
[K, f, loc, viol] = hifoomain(pars, order, init, obj, options);

% In the case of the H2 optimisation, the formulae used for norm and gradient computation is 
% valid only for the case where D22=0. If D22!=0 we have to apply a
% transformation on the controller to
% compensate for the use of the formulae for the shifted output y'= y-D22*u
% The field structure.D22 simply tells that the D22 term of the plants exists and
% the transformation should be applied. The transformation K->K1 is:
% K1.a = K.a - K.b*D22*(1+K.d*D22)^-1*K.c
% K1.b = K.b*(1 - D22*(1+K.d*D22)^-1*K.d)
% K1.c = (1 + K.d*D22)^-1*K.c
% K1.d = (1 + K.d*D22)^-1*K.d

%if we have H2 synthesis and D22 != 0 we change the controller 
if isstruct(options.structure) && isfield(options.structure,'D22')
    %check to see that D22 is not all zeros before making any changes
    if (sum(sum(options.structure.D22 ~= zeros(size(options.structure.D22)))) ~= 0)
        %we have to apply the transformation on the controller 
        N = K.d*options.structure.D22;
        I = eye(size(N));
        
        M = I + N;     % 1+Dk*D22
        cd = cond(M);      % the conditionment number for the patrix that should be inverted
        
        if cd*eps>1
            if cd == Inf    % 1+Dk*D22 == 0
                error('hifoo: the transformation matrix 1+Dk*D22 == 0, one of the closed-loop plants may not be well-posed \n');
            else
                warning('hifoo: the transformation matrix 1+Dk*D22 is singular to working precision, one of the closed-loop plants may not be well-posed \n');
            end
        end
        
        % if I arrive here it means that the matrix can be inverted
        
        K.d = M\K.d;
        K.c = M\K.c;
        
        K.a = K.a - K.b*options.structure.D22*K.c;
        K.b = K.b - K.b*options.structure.D22*K.d;
    end
end
        

if (ssinputplant || ssinputcont) % convert to ss object format of Control Systems Toolbox
    K = ss(K.a, K.b, K.c, K.d);
end
