function [f, g] = hifooobj(x, pars)

% THIS NEXT PART IS NOW DIFFERENT IN HIFOO 2.0
% MODIFIED by Marc Millstone
% if pars.objective = 'i', then each plant has its own objective in
% pars.plantinfo{i}.objective
% otherwise, pars.objective ='+' and we stabilize everything for the
% first part of hifoo, where we are searching for just a stabilizing
% controller
% Previously, in earlier versions, the hifoo optimization objective was
% %   pars.objective = 's': spectral abscissa
% %   pars.objective = 'r': inverse of complex stability radius
% %   pars.objective = 'p': inverse of complex stability radius
% %   pars.objective = 'h': H-infinity norm
% along with a penalty term on ||x||_2, controlled by pars.penalty,
% AND, in the case of the H-infinity norm, a barrier on
% the stability region, controlled by pars.barrier.
% The second is a special case of the third (with pars.barrier=0).
% The parameters define a reduced order model, with order pars.nhat.
% These parameters, normally known as Ahat, Bhat, Chat, Dhat, are encoded
% in a single parameter vector x for the purposes of optimization.
% The fixed data are defined in
%   pars.A, pars.B1, pars.B2, pars.C1, pars.C2, pars.D11, pars.D12, pars.D21
%
% It is not no longer assumed that D22 == 0.
%
% Modified by Marc Millstone to accept multiple plants and take the max
% over the norm of each plant.
% Pars now has a different structure as well.
% Modified: To allow norm constraints on plants.
% Modified to add strong stabilization support.  We must determine if
% the current plant is actually the controller.  IF so, we use the
% controller matrices, not the plant matrices.

% Modifed further by Georgia Deaconu to allow the H2 objective too.

% This function returns the optimization objective as well as the gradient
% with respect to the parameters. If the function turns out to be infinite,
% the gradient will be nan, which is correctly propogated by chainABCD.
% The SLICOT routine linorm is used if pars.slicot == 1; otherwise,
% the norm function in Matlab's control system toolbox is used.
% The function returns the norm of the controller as well.

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



mainobjective = pars.objective;
%First get which plants we are to optimize over
optimize = pars.optimize;
%This plants must meet a certain objective constraint
constraints = setdiff(1:numel(pars.plantinfo),optimize);
stabilize = pars.stabilize;

%We must get the function values for each plant in pars.plantinfo
fvals = zeros(1,numel(pars.plantinfo));  %To store function values
wvals = zeros(1,numel(pars.plantinfo)); % To store w vals for Hinf norm
slicot = pars.slicot;  % 1 if slicot is available, 0 otherwise
constraintPenalty = pars.penaltyConstraint;
f = 0;
g = 0;
nhat = pars.nhat;


structure = pars.structure;
% if there is a barrier term, we should save the w vals
% for it as well
if pars.barrier > 0
    wvals_barrier = zeros(1,numel(pars.plantinfo));
end
for i = 1:numel(pars.plantinfo)
    % We are stabilizing and this plant does not require
    % stabilization in this phase.
    if (mainobjective == '+' && isempty(find(stabilize == i)))
        continue;
    end
    
    
    
    plant = pars.plantinfo{i};
    
    if (isfield(plant,'controller') && plant.controller == 1)
        %This is for strong stabilization, this is the controller
        m = pars.m;
        p = pars.p;
        [Abig, Bbig, Cbig, Dbig] = getABCDhat(pars.nhat, m, p, x,pars.structure);
    else
        A = plant.A; B1 = plant.B1; B2 = plant.B2; C1 = plant.C1; C2 =plant.C2;
        D11 = plant.D11; D12 = plant.D12; D21 = plant.D21; D22 = plant.D22;
        structure = pars.structure;
        n = size(A,1);
        
        [Abig, Bbig, Cbig, Dbig] = ...
            getABCDbig(A, B1, B2, C1, C2, D11, D12, D21,D22,nhat, x,structure );
        
        m = pars.m;
        p = pars.p;
        [Ahat, Bhat, Chat, Dhat] = getABCDhat(nhat, m, p, x, structure);
    end
    %Check whether we are just stabilizing or allowing each plant to
    %have its own objective function
    if mainobjective == '+'
        objective = 's';
    else
        objective = plant.objective;
    end
    
    %if nargout < 2 % function value only,
    %commented out.  We need to get the function values anyway.
    if objective == 's' % spectral abscissa
        f = specabsc(Abig);
    elseif objective == 'r' % inverse of complex stability radius
        if (isfield(plant,'controller') && plant.controller == 1)
            I = eye(nhat);
            zero = zeros(nhat);
        else
            I = eye(n + nhat);
            zero = zeros(n + nhat);
        end
        [f,w] = hinfty(Abig, I, I, zero, slicot);
        wvals(i) = w;
    elseif objective == 'h' % most important case: H-infinity norm
        if isempty(Abig)
            % This is a SOF controller, so H-infinity norm is max(svd(Dbig))
            [S] = svd(Dbig);
            f = max(S);
            wvals(i) = 0;
        else
            [f,w] = hinfty(Abig, Bbig, Cbig, Dbig, slicot);
            wvals(i) = w;
        end
        
        % Now add the barrier term if necessary, ignore if SOF controller
        barrier = pars.barrier;
        if barrier > 0 && ~isempty(Abig)
            if (isfield(plant,'controller') && plant.controller == 1)
                I = eye(nhat);
                zero = zeros(nhat);
            else
                I = eye(n + nhat);
                zero = zeros(n + nhat);
            end
            % We will have to recompute this gradient below!
            [fbarrier,wbarrier] = hinfty(Abig, I, I, zero, slicot);
            f = f + fbarrier;
            wvals_barrier(i) = wbarrier;
            
        end
    elseif objective == 't'
        [f,gh2] = htwo(plant.A,plant.B1,plant.B2,plant.C1,plant.C2,plant.D11,...
            plant.D12,plant.D21,Ahat,Bhat,Chat,Dhat,pars.nvar,pars.structure);
    elseif objective == 'p'; % pseudospectral abscissa
        f = pseudospecabsc(Abig, pars.epsilon);
    else % cannot happen
        fprintf('hifooobj: invalid objective\n')
        f = nan;
    end
    
    %depending on whether the plant is part of the optimization or
    %constraint we set the fval differently
    
    if ((~isempty(find(optimize==i))))
        %Then the plant is in the optimize set
        fvals(i) = pars.weights(i)*f;
        
    elseif (~isempty(find(constraints == i)))
        %The plant is a constraint
        normconstraint = plant.constraint;
        fvals(i) = constraintPenalty*max(pars.weights(i)*f-normconstraint,0);
    end
end

if mainobjective == '+'
    %only care about stabilization
    [f, index] = max(fvals(stabilize));
    w = wvals(index);
else
    %Must account for constraint penalty
    [f, index] = max(fvals(optimize));
    if ~isempty(f)
        %We have to calcuate everything
        f = f+sum(fvals(constraints));
        w = wvals(optimize(index));
        if pars.barrier > 0
            wbarrier = wvals_barrier(optimize(index));
        end
    else
        %There are only constraints to be met
        f = sum(fvals(constraints));
    end
end


if nargout ==2
    % We compute the gradient of pars.plantinfo{optimize(index)}
    % in the case of the hinfinity norm, we use the f and w values
    % computed above
    
    g = 0;
    if mainobjective == '+'
        plants = stabilize(index);
    else
        plants = [optimize(index) constraints];
    end
    for i = 1:numel(plants)
        
        plant = pars.plantinfo{plants(i)};
        if (isfield(plant,'controller') && plant.controller == 1)
            %This is for strong stabilization
            
            m = pars.m;
            p = pars.p;
            [Abig, Bbig, Cbig, Dbig, nVarA, nVarB, nVarC,nVarD] = getABCDhat(pars.nhat, m, p, x,pars.structure);
        else
            
            A = plant.A; B1 = plant.B1; B2 = plant.B2; C1 = plant.C1; C2 =plant.C2;
            D11 = plant.D11; D12 = plant.D12; D21 = plant.D21; D22 = plant.D22;
            structure = pars.structure;
            n = size(A,1);
            nhat = pars.nhat;
            [Abig, Bbig, Cbig, Dbig] = ...
                getABCDbig(A, B1, B2, C1, C2, D11, D12, D21,D22, nhat, x,structure );
        end
        
        %Check whether we are just stabilizing or allowing each plant to
        %have its own objective function
        if mainobjective == '+'
            objective = 's';
        else
            objective = plant.objective;
        end
        % We have to calculate gradient if the plant corresponds to
        % maximum objective or when the constraint is not satisfied
        if i == 1 || ( i > 1 && fvals(plants(i)) > 0)
            if objective == 's'
                [f2, G] = specabsc(Abig);
                if (isfield(plant,'controller') && plant.controller == 1)
                    %This is for strong stabilization
                    if ~isstruct(structure) && structure == -1
                        gtemp = reshape(real(G),nhat*nhat,1);
                        gtemp = pars.strongstabweight*[gtemp; zeros(pars.nvar-nhat*nhat,1)];
                    else
                        gtemp = G(structure.Ahat);
                        gtemp = pars.strongstabweight*[gtemp; zeros(nVarB+nVarC+nVarD,1)];
                    end
                else
                    gtemp = chainABCD(n, nhat, B2, C2, D12, D21, D22,G, ...
                        zeros(size(Bbig)), zeros(size(Cbig)), zeros(size(Dbig)),x,structure);
                end
            elseif objective == 'r'
                if (isfield(plant,'controller') && plant.controller == 1)
                    I = eye(nhat);
                    zero = zeros(nhat);
                    zb = zeros(nVarB,1);
                    zc = zeros(nVarC,1);
                    zd = zeros(nVarD,1);
                else
                    I = eye(n + nhat);
                    zero = zeros(n + nhat);
                end
                [f2,w,G] = hinfty(Abig, I, I, zero, slicot,f,w);
                if (isfield(plant,'controller') && plant.controller == 1)
                    if ~isstruct(structure) && structure == -1
                        gtemp = pars.strongstabweight*[reshape(real(G),numel(G),1);...
                            zb(:); zc(:); zd(:)];
                    else
                        gtemp = real(G(structure.Ahat));
                        gtemp = pars.strongstabweight*[gtemp;...
                            zb(:); zc(:); zd(:)];
                    end
                else
                    gtemp = chainABCD(n, nhat, B2, C2, D12, D21,D22, G, ...
                        zeros(size(Bbig)), zeros(size(Cbig)), zeros(size(Dbig)),x,structure);
                end
                
            elseif objective == 'p'
                [f2, G] = pseudospecabsc(Abig, pars.epsilon);
                gtemp = chainABCD(n, nhat, B2, C2, D12, D21, D22,G, ...
                    zeros(size(Bbig)), zeros(size(Cbig)), zeros(size(Dbig)),x,structure);
            elseif objective == 'h'
                if isempty(Abig)
                    %SPecial case:  This is a SOF controller
                    % Just compute the gradient here
                    [U,S,V] = svd(Dbig);
                    u = U(:,1);  % left singular vector for largest singular value
                    v = V(:,1);  % right singular vector for largest singular value
                    uvt = u*v';
                    if ~isstruct(structure) && structure == -1
                        gtemp = pars.strongstabweight*[reshape(real(uvt),numel(uvt),1)];
                    else
                        uvt = uvt(structure.Dhat);
                        gtemp = pars.strongstabweight*[reshape(real(uvt),numel(uvt),1)];
                    end
                else
                    [f2,w, GA, GB, GC, GD] = hinfty(Abig, Bbig, Cbig, Dbig, slicot,f,w);
                    
                    if (isfield(plant,'controller') && plant.controller == 1)
                        %This is for strong stabilization
                        if ~isstruct(structure) && structure == -1
                            gtemp = pars.strongstabweight*[reshape(real(GA),numel(GA),1); reshape(real(GB),numel(GB),1);reshape(real(GC),numel(GC),1);...
                                reshape(real(GD),numel(GD),1)];
                        else
                            GA = GA(structure.Ahat);
                            GB = GB(structure.Bhat);
                            GC = GC(structure.Chat);
                            GD = GD(structure.Dhat);
                            gtemp = pars.strongstabweight*[reshape(real(GA),numel(GA),1); reshape(real(GB),numel(GB),1);reshape(real(GC),numel(GC),1);...
                                reshape(real(GD),numel(GD),1)];
                        end
                    else
                        gtemp = chainABCD(n, nhat, B2, C2, D12, D21,D22, GA, GB, GC, GD,x, structure);
                    end
                    
                    %now compute the gradient from the barrier term if utilized
                    if pars.barrier > 0
                        if (isfield(plant,'controller') && plant.controller == 1)
                            I = eye(nhat);
                            zero = zeros(nhat);
                            zb = zeros(nVarB,1);
                            zc = zeros(nVarC,1);
                            zd = zeros(nVarD,1);
                        else
                            I = eye(n + nhat);
                            zero = zeros(n + nhat);
                        end
                        [f2,w,G] = hinfty(Abig, I, I, zero, slicot,f,wbarrier);
                        if (isfield(plant,'controller') && plant.controller == 1)
                            if ~isstruct(structure) && structure == -1
                                gtemp = pars.strongstabweight*[reshape(real(G),numel(G),1);...
                                    zb(:); zc(:); zd(:)];
                            else
                                gtemp = real(G(structure.Ahat));
                                gtemp = pars.strongstabweight*[gtemp;...
                                    zb(:); zc(:); zd(:)];
                            end
                        else
                            gtemp_bar = chainABCD(n, nhat, B2, C2, D12, D21,D22, G, ...
                                zeros(size(Bbig)), zeros(size(Cbig)), zeros(size(Dbig)),x,structure);
                        end
                        gtemp = gtemp + gtemp_bar;
                    end
                    
                end
            elseif objective == 't'
                gtemp = gh2;
            else % cannot happen
                fprintf('hifooobj: invalid objective\n')
                f = nan;
                gtemp = nan*ones(pars.nvar,1);
            end
            
            
            
            
            %The optimize plant ALWAYS comes first, so we can set it the first time.
            if (i == 1|| mainobjective == '+')
                %Then the plant is in the optimize set or we are just stabilizing
                g = pars.weights(plants(i))*gtemp;
            else
                %The plant is a constraint
                g = g+pars.weights(plants(i))*constraintPenalty*gtemp;
            end
        end
    end
end





% the penalty term is to prevent the norm of the controller from blowing up
penalty = pars.penalty;
if penalty > 0
    xnorm = norm(x);
    if xnorm > 0
        f = f + penalty*xnorm;
        if nargout >= 2
            g = g + penalty*x/xnorm; % gradient of the 2-norm
        end
    else
        d = randn(pars.nvar);  % random subgradient in this case
        if nargout >= 2
            g = g + penalty*d/dnorm;
        end
    end
end
% the barrier term addresses the fact that the H-infinity norm does
% not necessarily increase, and in fact may decrease, as the eigenvalues
% of Abig go to zero, unlike the inverse of the stability radius.
% This makes the objective very difficult to optimize, as the value
% at the stability boundary is by definition infinity.  Adding a small
% multiple of the inverse of the stability radius to the H-infinity
% objective alleviates this difficulty, as it blows up as the
% stability boundary is approached, behaving as a barrier to instability.

% THIS CODE HAS BEEN MOVED ABOVE, SO WE CAN HANDLE MULTIPLE PLANTS

% objective = pars.objective;
% if objective == 'h'
%     barrier = pars.barrier;
%     if barrier > 0
%         I = eye(n + nhat);
%         if nargout < 2
%             stabradinv = hinfty(Abig, I, I, zeros(n + nhat), slicot);
%             f = f + barrier*stabradinv;
%         else
%             I = eye(n + nhat);
%             [stabradinv, G] = hinfty(Abig, I, I, zeros(n + nhat), slicot);
%             f = f + barrier*stabradinv;
%             g = g + barrier*chainABCD(n, nhat, B2, C2, D12, D21,D22, G, ...
%               zeros(size(Bbig)), zeros(size(Cbig)), zeros(size(Dbig)),x,structure);
%         end
%     end
% end
