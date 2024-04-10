function [pars, options] = processOptions(pars, options)
% Process the option input 
% This will also convert new option names to the original, internal Hifoo 
% names

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

order = pars.nhat;
m = pars.m;
p = pars.p;
pars.optimize = options.optimize;
pars.stabilize = options.stabilize;
pars.penaltyConstraint = 0;



%FIRST TRANSLATE OPTIONS FROM NEW NAMES TO OLD NAMES

%  weighNormK --> penalty
%  augmentHinf --> barrier
%  rho.init --> penaltyConstraintInit
%  rho.max  --> penaltyConstraintMax
%  rho.multiply --> penaltyConstraintMultiply

if isfield(options,'weightNormK')
    options.penalty = options.weightNormK;
    options = rmfield(options,'weightNormK');
end
if isfield(options,'augmentHinf')
    options.barrier = options.augmentHinf;
    options = rmfield(options,'augmentHinf');
end



if isfield(options, 'rho')
    rho = options.rho;
    options.penaltyConstraintInit = rho.init;
    options.penaltyConstraintMultiply = rho.multiply;
    options.penaltyConstraintMax = rho.max;
    options = rmfield(options, 'rho');
end



nvar = (order + m)*(order + p);

% Structure is alredy set up and verified in hifoo.
if isstruct(options.structure)
    %We only optimize over a subset of controller entries
    varsD = m*p;
    
    %if we have one matrix we have them all, as set in processStructure
    if (isfield(options.structure,'Ahat') && isfield(options.structure,'Bhat')...
            && isfield(options.structure,'Chat') && isfield(options.structure,'Dhat'))
        
        varsA = size(options.structure.Ahat,1);
        varsB = size(options.structure.Bhat,1);
        varsC = size(options.structure.Chat,1);
        varsD = size(options.structure.Dhat,1);
        
        if varsD == 0       %this means that we have unique solution for Dk
            varsD = m*p;            %this is ugly but works.   
        end
        nvar = varsA + varsB + varsC + varsD;
    end
    
    if isfield(options.structure,'Dk') || isfield(options.structure,'V')

        nvar = nvar - varsD;
        
        if isempty(options.structure.Dk)    %it means we have stuff to optimize
            nvar = nvar + size(options.structure.V,2);
        end
    end
end

pars.nvar = nvar;
pars.structure = options.structure;  %already verified and set up in hifoo.m
% hifoo options (no call to setdefaults)

 
if ~isfield(options,'fast')
    options.fast = 0;
else
    if options.fast ~= 0 && options.fast ~= 1
        error('hifoo: options.fast must be either 0 or 1')
    end
end
    
if ~isfield(options, 'penalty') 
    options.penalty = 0;
elseif isempty(options.penalty)
    options.penalty = 0;
elseif ~(isposreal(options.penalty) | options.penalty == 0)
    error('hifoo: input "options.penalty" must be a nonnegative real scalar')
end
penalty = options.penalty; % penalty is copied to pars.penalty below


if ~isfield(options,'weights') 
    pars.weights = ones(1,pars.nplants);
else
    % verify that weights are strictly positive
    if length(options.weights) ~= pars.nplants 
        error('hifoo: input options.weight is of the wrong dimension\n')
    elseif sum(options.weights <=0) ~= 0
        error('Hifoo:  Weights must be strictly positive')
    else
        pars.weights = options.weights;
    end
end
    


if ~isfield(options, 'barrier') 
    options.barrier = 0;
elseif isempty(options.barrier)
    options.barrier = 0;
elseif ~(isposreal(options.barrier) | options.barrier == 0)
    error('hifoo: input "options.barrier" must be a nonnegative real scalar')
end
if ~isfield(options, 'cpumax') % quit when cpu time in seconds exceeds this
    options.cpumax = inf;
elseif ~isposreal(options.cpumax)
    error('hifoo: input "options.cpumax" must be a positive scalar')
end
if ~isfield(options, 'normtol')
    options.normtol = 1e-3; % larger default than HANSO default
elseif ~isposreal(options.normtol)
    error('hifoo: input "options.normtol" must be a positive scalar')
end 
if ~isfield(options, 'evaldist') 
    options.evaldist = 1e-3; 
elseif ~isposreal(options.evaldist)
    error('hifoo: input "options.evaldist" must be a positive scalar')
end
if ~isfield(options, 'nrand') % number of starting points besides init
    options.nrand = 3;
elseif ~isposint(options.nrand)
    % nrand is 0 is not allowed: if init came from output of a hifoo run,
    % perhaps a different order, objective is likely not to be smooth
    % there and BFGS may fail immediately: some perturbation is needed
    error('hifoo: input "options.nrand" must be a positive integer')
end
if ~isfield(options,'penaltyConstraintInit')
    options.penaltyConstraintInit = 100;
end
if ~isfield(options,'penaltyConstraintMultiply')
    options.penaltyConstraintMultiply = 100;
end

if ~isfield(options,'penaltyConstraintMax')
    options.penaltyConstraintMax = 1e10;
end

if ~isfield(options,'hinfalg')
    options.hinfalg = 'matlab';
elseif ~strcmp(options.hinfalg,'matlab') && ~strcmp(options.hinfalg,'slicot')
    error('hifoo: options.hinfalg must be either matlab or slicot')
end

if ~isfield(options,'hinftol')
    options.hinftol = 1e-7;
else ~isposreal(options.hinftol)
    error('hifoo: options.hinftol must be a positive scalar ')
end


if ~isfield(options,'maxit')
    options.maxit = 1000;
elseif ~isposint(options.maxit)
    error('hifoo: options.maxit must be a positive integer')
end     

pars.barrier = options.barrier;
