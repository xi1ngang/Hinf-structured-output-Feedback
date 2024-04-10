function [K, f, loc, viol] = hifoomain(pars, order, init, obj, options)
% main routine for HIFOO: H-infinity fixed-order optimization
% intended to be called by hifoo only
% see comments in hifoo.m
% pars and init use the structure format (other formats already converted)
%
%
% The output gives the controller, the objective function, a certificate
% as to the quality of the minimizer and a vector of constraint violations



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

%All plants must have same dimensions here.  So just grab the first plant
%It has already been checked that the dimensions of all subsequent plants
%match
m = pars.m; % number of cols of B2 (number of rows of Chat)
p = pars.p; % number of rows of C2 (number of cols of Bhat)
pars.nhat = order;

% function that computes H-infinity norm, H2 norm, complex stability
% radius,
% pseudospectral abscissa or spectral abscissa, depending on pars.obj
pars.fgname = 'hifooobj';

[pars,options] = processOptions(pars,options);
viol = Inf*ones(1,numel(setdiff(1:numel(pars.plantinfo),pars.optimize)));
cpufinish = cputime + options.cpumax;
penalty = options.penalty;  %Moved into pars below
barrier = options.barrier; % barrier is copied to pars.barrier below
nrand = options.nrand;
prtlevel = options.prtlevel; % set in hifoo
if prtlevel > 0 % intro line already printed by hifoo
    disp(options)
end

%Not important for nultiple plants
% if pars.D22 == 0
%     fprintf('hifoo: D22 = 0, exactly.  Using simplified closed loop equations\n')
% else
%     fprintf('hifoo: D22 not exactly zero, using nonlinear closed loop equations\n')
% end


if prtlevel > 0
    fprintf('hifoo: number of design variables is %d\n', pars.nvar)
end
options.prtlevel = max(0, prtlevel-1); % to reduce output from HANSO



% for pseudospectral abscissa, copy epsilon info from options to pars
if ~isempty(find(obj == 'p'))
    if isfield(options,'epsilon') % defines epsilon-pseudospectral abscissa
        if ~isposreal(options.epsilon)
            error('hifoo: input "options.epsilon" must be a positive real scalar')
        else
            pars.epsilon = options.epsilon; %  default value, for passing to hifooobj via PARS, not options
        end
    else
        pars.epsilon = .01;
    end
    if exist('pspa') ~= 2
        error('hifoo: to compute pseudospectral abscissa, "pspa" must be installed from www.cs.nyu.edu/overton/software/pspa/')
    end
end
% check out what software is available to compute H-infinity norm
if ~isempty(find(obj == 'h')) | ~isempty(find(obj == 'r'))
    if strcmp(options.hinfalg,'slicot') && exist('linorm') == 3 % SLICOT is the fastest choice (3 means mex file)
        slicot = 1;
        if prtlevel > 0
            fprintf('hifoo: using "linorm" from SLICOT to compute H-infinity norm and complex stability radius\n')
        end
    elseif exist('lti') == 2% Control System Toolbox is the default (2 means m-file)
        if prtlevel > 0
            fprintf('hifoo: using Control System Toolbox to compute H-infinity norm and complex stability radius')
        end
        slicot = 0;  % Control System Toolbox
    else
        error('hifoo: either Control System Toolbox  or SLICOT (linorm) must be installed when input "obj" is "h" or "r"')
    end
else
    slicot = 0; % don't need H-infinity norm or complex stability radius
end
pars.slicot = slicot;

% check out what software is available to compute H2 norm
if ~isempty(find(obj == 't'))
    if exist('lyap') ~= 2
        error('hifoo: lyap function from the Control System Toolbox must be installed when input "obj" is "t"');
    end
end

% check whether HANSO is available: must be version 2.0 or higher which does not need quadprog   
if exist('hanso') == 2
    if prtlevel > 0
        fprintf('hifoo: using "hanso" for optimization (should be version 2.0 or higher)\n')
    end
else
    error('hifoo: "hanso (version 2.0)" must be installed from www.cs.nyu.edu/overton/software/hanso/hanso2_0')
end
% initialization
if isempty(init)
    % no initial point provided, starting points are randomly generated
    x0 = randn(pars.nvar,nrand);
    if prtlevel > 0
        if nrand > 1
            fprintf('hifoo: no initial point provided, using %d randomly generated starting points\n', nrand)
        else
            fprintf('hifoo: no initial point provided, using a randomly generated starting point\n')
        end
    end
else
    structure = options.structure;
    % set up initial point for optimization, expanding (or truncating) the
    % given data to correspond to the desired order if necessary
    xinit = xinitsetup(m, p, order, init, obj, pars,structure);
    
    % starting points include progressively bigger perturbations of xinit
    scalepert = norm(xinit)*logspace(-3, -1, nrand); % spaced from 0.001 to 0.1
    xinitpert = randn(pars.nvar, nrand)*diag(scalepert);
    x0 = [xinit    xinit*ones(1,nrand) + xinitpert];
    if prtlevel > 0
        if nrand > 1
            fprintf('hifoo: supplementing initial point with %d randomly generated perturbations\n', nrand)
        else
            fprintf('hifoo: supplementing initial point with a randomly generated perturbation\n')
        end
    end
end
oldopts = options;
% in hifoo 1.0, obj was just one character
% in hifoo 2.0 or higher, it might be just one character, or it might be a string
% if it is one character, one of the following may be "true", but we could
% just as well say "if ~isempty(options.stabilize)", since this will also
% be "true" in that case.  Therefore the following line, used in hifoo 2.0,
% is commented out and replaced by the simpler test.


 % if obj == 'h' | obj == 'r' | obj == '+' | ~isempty(options.stabilize)


if ~isempty(options.stabilize)
    % for these objectives, first step is to get stable
    % If there is more than one plant, we always stabilize all plants first
    pars.objective = '+';   % to this end, minimize spectral abscissa
    pars.penalty = 0;       % without any penalty term (barrier irrelevant)
    options.fvalquit = 0;   % stopping when get stable
    options.maxit = 1000;   % run a long time if necessary
    if prtlevel > 0
        fprintf('hifoo: searching for stabilizing order %d controllers\n', order)
        if nrand > 1
            fprintf('hifoo: if stabilization time is excessive, reduce options.nrand\n')
        end
    end
    % use BFGS alone at first, as in most cases this will be sufficient
    % call it one starting point at a time, because when objective is '+',
    % want to stop when a stable point is found, while for objectives
    % 'h', 't' and 'r', want to find several stable points for initializing 
    % optimization
    pars_stabilize = pars;
    pars_stabilize.optimize = options.stabilize;
    for k=1:size(x0,2)
        options.x0 = x0(:,k);
        % might waste too much of the allotted cputime stabilizing many
        % start points but complicated to avoid, user can adjust options
        options.cpumax = cpufinish - cputime;
        [xB(:,k),fB(k),gB(:,k)] = bfgs(pars_stabilize,options);
        if cputime > cpufinish
            if prtlevel > 0
                fprintf('hifoo: quit stabilizing since CPU time limit exceeded\n')
            end
            break % not return
        end
        if obj == '+' & fB(k) < 0 % only want one stable point
            break % not return
        end
    end
    % if BFGS did not find a stable point, try gradient sampling,
    % initializing with the best half of the points found by BFGS,
    % starting with the lowest
    % (no point using local bundle, not trying to verify optimality)
    if cputime < cpufinish & min(fB) >= 0
        if options.fast == 1
            if prtlevel > 1
                fprintf('hifoo: skipping gradient sampling\n');
            end
        else
            if prtlevel > 1
                fprintf('hifoo: BFGS did not find any stabilizing controllers, trying gradient sampling\n')
            end
            [fBsort,indx] = sort(fB);
            indx = indx(1:ceil(length(fB)/2));
            xBsort = xB(:,indx);
            options.maxit = 100; % gradient sampling is more expensive than BFGS
            for k = 1:length(indx)
                options.x0 = xBsort(:,k);
                options.cpumax = cpufinish - cputime;
                [xGS, fGS, gGS] = gradsamp(pars_stabilize,options);
                if cputime > cpufinish
                    if prtlevel > 0
                        fprintf('hifoo: quit stabilizing since CPU time limit exceeded\n')
                    end
                    break % not return
                end
                if fGS < 0 % settle for one stable point
                    xB = xGS; fB = fGS; gB = gGS;
                    break % not return
                end
            end
        end
    end
    stabindx = find(fB < 0);
    nstabpts = length(stabindx);
    if nstabpts > 0
        if prtlevel > 0
            if min(fB) > -1e-8
                qualification = '(barely) ';
            else
                qualification = '';
            end
            if nstabpts == 1
                fprintf('hifoo: found a %sstabilizing controller', qualification)
            else
                fprintf('hifoo: found %d %sstabilizing controllers', nstabpts, qualification)
            end
            if obj ~= '+'
                fprintf(' for initializing optimization\n')
            else
                fprintf(' , quitting\n')
            end
        end
        x0 = xB(:,stabindx);  % at least one stable point was found
        f0 = fB(:,stabindx);
        foundstablepoint = 1;
    else
        if prtlevel > 0
            fprintf('hifoo: could not find a stabilizing order %d controller\n', order)
            fprintf('hifoo: returning controller with best spectral abscissa %g instead\n', min(fB))
            fprintf('hifoo: try specifying a new "init" or increasing "options.cpumax"\n')
        end
        foundstablepoint = 0;
    end
    if obj == '+'| ~foundstablepoint % nothing else to do
        % fB and xB, not f0 and x0, as f0 and x0 are empty in 2nd case
        [f, k] = min(fB); % in both cases, return one with lowest abscissa
        if obj ~= '+'
            f = inf; % but return infinite objective when obj is 'h', 't' or 'r'
        end
        x = xB(:,k);
        K = [];
        [K.a, K.b, K.c, K.d] = getABCDhat(order, m, p, x, options.structure);
        loc.dnorm = nan;
        loc.evaldist = nan;
        return
    end
end
options=oldopts;

% i here means independent or individual.  hifooobj will look at the set objective in pars.plantinfo
pars.objective = 'i'; % No one objective function, each plant carries its own
objectivename = objname(obj, pars);
pars.objname = objectivename; % may be handy for display purposes
pars.penalty = penalty;   % penalty term on ||x||_2
if penalty > 0
    penaltystring = ' (plus penalty term)';
else
    penaltystring = '';
end
pars.barrier = barrier;   % barrier term: multiple of inverse of stability radius
if barrier > 0
    barrierstring = ' (plus barrier term)';
else
    barrierstring = '';
end
options.fvalquit = -inf; % do not quit early
options.x0 = x0;   % multiple starting points, stable if obj is 'h', 't' or 'r'
if prtlevel > 0
    fprintf('hifoo: optimizing %s%s%s for order %d controller\n', ...
        objectivename, penaltystring, barrierstring, order)
    fprintf('hifoo: if optimization time is excessive, reduce options.cpumax and/or options.nrand\n')
end
% run HANSO from all the valid starting points, returning only best result
% repeat runnign hanso with larger and larger penalty mutlipliers
% until we find a solution which satisfies all constraints
CONSTRAINTS_SATISFIED = 0;
pars.penaltyConstraintInit = options.penaltyConstraintInit;
pars.penaltyConstraint = options.penaltyConstraintInit;
while ~CONSTRAINTS_SATISFIED && pars.penaltyConstraint <= options.penaltyConstraintMax
    options.cpumax = cpufinish - cputime; % time left after stabilization
    cpufinish = cputime + options.cpumax;
    % in HIFOO 2.0, called BFGS (part of HANSO 1.0) directly 
    % in HIFOO 3.5, prefer to call HANSO 2.0, as now even the BFGS phase
    % alone returns a meaningful LOC in the nonsmooth case: no longer need
    % to use local bundle or gradient sampling
    % Note: HANSO, unlike BFGS, returns only the best result, not one for each
    % starting point
    options.samprad = []; % turn off gradient sampling (HANSO 2.0)
    % in case of HANSO 1.0 is called, turn off local bundle and gradient sampling
    options.phasenum = [size(options.x0,2) 0 0]; 
    [xsol,f,loc] = hanso(pars,options);
    

    
    %Now we check constraints to see if they are satisfied
    constraints = setdiff(1:numel(pars.plantinfo),pars.optimize);

	%       We need to compute *only* the value of the constraint
	% 		infeasibilities. We do this by creating a new parameter collection
	%		and set the set of plants to be only the set of plants that appear 
	%		in the constraint list.  The set of plants to optimize is empty.
	
    parsToEvaluateOnlyConstraints = pars;
    parsToEvaluateOnlyConstraints.plantinfo = {pars.plantinfo{constraints}};
    parsToEvaluateOnlyConstraints.optimize  = [];

	% This will evaluate only the constraint violations
    fconstraint = hifooobj(xsol,parsToEvaluateOnlyConstraints);
    if fconstraint == 0 && options.fast == 0 
        if prtlevel > 1
            fprintf('hifoo: attempting to improve solution via gradient sampling\n');
        end
        % RUN GRADIENT SAMPLING (via HANSO) TO IMPROVE SOLUTION FROM BFGS ABOVE
        % Only run from the best point above.
        options.x0 = xsol;
        options.cpumax = cpufinish - cputime;
        % options.phasenum = [0,1,2]; was used by HANSO 1.0 to tell it to skip BFGS.
        % not used by HANSO 2.0, which will try BFGS and then immediately
        % switch to gradient sampling when BFGS fails as it has already been run
        % But gradient sampling is not done if the number of variables > 100 
        % as it is too expensive: then nothing will happen so xsol will not change
        options = rmfield(options, 'samprad'); % turn gradient sampling back on
        options = rmfield(options,'phasenum');  % turn gradient sampling back on (hanso 1.0)
        [xsol, f, loc] = hanso(pars,options);
    end
    
    % We must verify that constraints  are still
    % satisfied, as solution may not longer be feasible
    % after gradient sampling.
    fconstraint = hifooobj(xsol,parsToEvaluateOnlyConstraints);
    


    if fconstraint  == 0
        %CONSTRAINTS SATISFIED
        CONSTRAINTS_SATISFIED = 1;
    else
        CONSTRAINTS_SATISFIED = 0;
        %constraints not satisfied yet, increase penalty
        if isempty(pars.optimize)
            % increasing penalty will not effect the outcome, constraint cannot
            % be satisfied
            if prtlevel > 0
                fprintf('hifoomain: Optimize set is empty and constraints are not satisfied. Quitting...\n')
            end
            break;
        else
            
            pars.penaltyConstraint = options.penaltyConstraintMultiply*pars.penaltyConstraint;
            
            if prtlevel > 0
                fprintf('hifoomain: Constraints were not satisfied, increasing constraint penalty to %d\n',...x
                    pars.penaltyConstraint);
            end
            % Set initial points to be final points of previous phase
            %                 options.x0 = xsol;
            
        end
    end
    
    
    
    if cputime > cpufinish
        if prtlevel > 0
            fprintf('hifoomain: cpu time exceeded.  Exiting. Constraints may not be satisfied\n')
        end
        break;
    end
end

% Compute the constraint violations of each term in the constraint set
constraints = setdiff(1:numel(pars.plantinfo),pars.optimize);
viol = zeros(1,numel(constraints));

% Here yet again, we have the issue that originally I combined
% 	   the constraint computations with the hifooobj.  We need to 
%  	    compute the violation of each constraint, c_i>0.  To do
% 	    this we create a new parameter list, parsToEvaluateOnlyConstraints with only 1 plant.
%	    We then loop over each plant in the constraint set to compute
%	    the needed violation.  
for i = 1:numel(constraints)
    parsToEvaluateOnlyConstraints = pars;
    parsToEvaluateOnlyConstraints.nplants = 1;
    parsToEvaluateOnlyConstraints.plantinfo = {pars.plantinfo{constraints(i)}};
    parsToEvaluateOnlyConstraints.penaltyConstraint = 1;
    parsToEvaluateOnlyConstraints.optimize = [];
    viol(i) = hifooobj(xsol,parsToEvaluateOnlyConstraints);
end


% The final objective function above includes any terms from violated constraints.
% We recompute the final objective value with only terms in the actual objective portion
% of the problem.

parsToEvaluateObjective = pars;
parsToEvaluateObjective.nplants = length(pars.optimize);
parsToEvaluateObjective.plantinfo = {pars.plantinfo{pars.optimize}};

f = hifooobj(xsol,parsToEvaluateObjective);

% Verify and issue a warning if any of the closed loop systems are
% approaching the stability boundary
verifyStabilityBoundary(xsol, pars,prtlevel);

if isnan(loc.dnorm)
    loc.dnorm = zeros(size(loc.dnorm));
end

if prtlevel > 0
    if ~isempty(pars.optimize)
        
        fprintf('hifoo: best order %d controller found has %s%s%s %g\n', ...
            order, objectivename, penaltystring, barrierstring, f)
    else
        fprintf('hifoo: returning order %d controller attempting to satisfy given contraints\n', order)
    end
    if exist('loc')
        fprintf(' with local optimality measure: dnorm = %5.1e, evaldist = %5.1e\n',...
            loc.dnorm, loc.evaldist')
    end
end

% return the best point found, translated into controller format
K = [];
[K.a, K.b, K.c, K.d] = getABCDhat(order, m, p, xsol,pars.structure);
