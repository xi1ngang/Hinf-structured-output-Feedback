function [structure_out] = processStructure(pars,options, order)
% Verify that the given options.structure is of the correct format
% Convert from User format to internal format for Hifoo

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

p = pars.p;
m = pars.m;

%changed by Georgia Deaconu -> the structure is converted into a matrix that
%allows to extract only the elements that should be non-zero in the
%controller

% Convert between new and old naming scheme
if isfield(options,'struct')
    options.structure = options.struct;
    options = rmfield(options,'struct');
    structure = options.structure;
    %structure exists, convert from .a, .b, .c, .d to .Ahat, .Bhat, .Chat, .Dhat
    if isfield(structure,'a')
        structure.Ahat = structure.a;
        structure = rmfield(structure,'a');
    end
    if isfield(structure,'b')
        structure.Bhat = structure.b;
        structure = rmfield(structure,'b');
    end
    if isfield(structure,'c')
        structure.Chat = structure.c;
        structure = rmfield(structure,'c');
    end
    if isfield(structure,'d')
        structure.Dhat = structure.d;
        structure = rmfield(structure,'d');
    end
    options.structure = structure;
end


if ~isfield(options,'structure')
   %For internal use, no specified structure so set this to -1
   structure = -1;
else
    % The user has specified a structure so we must check its validity.
    % check for the old type of structure
    r = []; oldStruct = 0;
    if isfield(structure,'Ahat')
        r = [r; structure.Ahat(:)];
    end
    if isfield(structure,'Bhat')
        r = [r; structure.Bhat(:)];
    end
    if isfield(structure,'Chat')
        r = [r; structure.Chat(:)];
    end
    if isfield(structure,'Dhat')
        r = [r; structure.Dhat(:)];
    end
    
    if max(isnan(r)) < 1 % no NaN in any of the structure matrix so we are dealing with the old type of structure
        oldStruct = 1;
    end	
    
    % The dimension checks for each matrix were left unchanged, the new part
    % is the call to the function getStructureMatrix 
    % [structure.Ahat, structure.Afix] = getStructureMatrix(structure.Ahat,oldStruct) results in a matrix Ahat
    % that no longer has the same dimensions as the controller but that allows us to 
    % extract only the elements that should be optimized : xA = sturcture.Ahat*K.a extracts 
    % the elements to be opimized from the A matrix of the controller
    % Afix contains the values of the parameters that shouldn't be changed
    
    structure = options.structure;
    %AHAT  
    if ~isfield(structure,'Ahat')
        structure.Ahat = eye(order*order);      %all the elements should be nonzero so the selection matrix is I 
	structure.Afix = zeros(order,order);    % the fixed part is zero
        n1 = order; n2 = order;
    else 
        %Check dimensions
        [n1 n2] = size(structure.Ahat);
        if n1 ~= order
               error('Error in Structure: Ahat dimension does not match order');   
        elseif n1 ~= n2
            error('Error in structure: Ahat is not square');
        end
        
	%Check if we have Infs in the structure
        if max(max(isinf(structure.Ahat))) > 0
            error('Error in structure: Ahat cannot contain Inf elements');
        end

        %Check type, make sure we only have ones and zeros
        if ~islogical(structure.Ahat)
            if(oldStruct && (max(max(structure.Ahat)) > 1) | (min(min(structure.Ahat)) < 0) )
                error('Error in structure: using old structure Ahat cannot be converted to type logical as it contains non-binary elements')
            end
        end
        
        %Perform the transformations on the structure matrix
        [structure.Ahat,structure.Afix] = getStructureMatrix(structure.Ahat,oldStruct);
    end
    
    %BHAT
     if ~isfield(structure,'Bhat')
        structure.Bhat = eye(order*p);
        n3 = order; p1 = p;
	structure.Bfix = zeros(n3,p1);
    else 
        %Check dimensions
        [n3 p1] = size(structure.Bhat);
      
        if n1 ~= n3
            error('Error in structure: Ahat and Bhat do not have compatible dimensions');
        end
        %Check type
        if ~islogical(structure.Bhat)
            if ( oldStruct && (max(max(structure.Bhat)) > 1) | (min(min(structure.Bhat)) < 0))
                error('Error in structure: using old structure Bhat cannot be converted to type logical as it contains non-binary elements')
            end
        end
        
        [structure.Bhat, structure.Bfix] = getStructureMatrix(structure.Bhat,oldStruct);
     end
    
      %CHAT   
      if ~isfield(structure,'Chat')
          structure.Chat = eye(m*order);
          m1 = m; n4 = order;
	  structure.Cfix = zeros(m1,n4);
      else
          [m1 n4] = size(structure.Chat);
           if  n1 ~= n4
            error('Error in structure: Ahat and Chat do not have compatible dimensions');
           end
        if ~islogical(structure.Chat)
            if ( oldStruct && (max(max(structure.Chat)) > 1) | (min(min(structure.Chat)) < 0))
                error('Error in structure: using old structure Chat cannot be converted to type logical as it contains non-binary elements')
            end
        end
        
        [structure.Chat,structure.Cfix] = getStructureMatrix(structure.Chat,oldStruct);
      end
      
      %DHAT 
      if ~isfield(structure,'Dhat')
          structure.Dhat = eye(m*p);
          Z = [];
          m2 = m; p2 = p; 
	  structure.Dfix = zeros(m2,p2);
      else
         [m2 p2] = size(structure.Dhat);
         if  any(size(structure.Dhat) ~= [m p])
             error('Error in structure: Dhat has incorrect dimensions');
         elseif n1 ~= n4
            error('Error in structure: Ahat and Chat do not have compatible dimensions');
         elseif p1 ~= p2
            error('Error in structure: Bhat and Dhat do not have compatible dimensions');
         elseif m1~=m2 
            error('Error in structure: Chat and Dhat do not have compatible dimensions');
         end
         
         if ~islogical(structure.Dhat)
            if (oldStruct && (max(max(structure.Dhat)) > 1) | (min(min(structure.Dhat)) < 0))
                error('Error in structure: using old structure Dhat cannot be converted to type logical as it contains non-binary elements')
            end
         end

         [structure.Dhat,structure.Dfix,Z] = getStructureMatrix(structure.Dhat,oldStruct);
      end
            
end

%construct the big system that gives the restrictions on Dhat
% For multiple plants H2 optimization the constraint 
% D11 + D12*K.d*D21 = 0 must be verified for all the plants
% we construct the equivalent system where D11 holds the D11 matrix from
% all the plants and so on

j = 1;
for i = 1:numel(pars.plantinfo)
	thisplant = pars.plantinfo{i};

	if (thisplant.objective == 't')
		m = size(thisplant.B2,2);
		p = size(thisplant.C2,1);
		D11{j} = thisplant.D11;
		D12{j} = thisplant.D12;
        D21{j} = thisplant.D21;
        D22{j} = thisplant.D22;
		j = j+1;
	end
end

if (j > 1)
    % If we are in the case of H2 opimization (multiple plants eventually) we have
    % to make sure that all the plants have equal D22. If the D22 matrix
    % are different, the transformation that we apply at the end on the
    % controller is not valid (in hifoo.m)

    D = D22{1};
    for i = 2:j-1
        if sum(sum(D ~= D22{i})) ~= 0
            error('Error: the plants with H2 objective cannot have different D22 terms \n');
        end
    end 
    
    A = [];
    b = [];
    
    % Construct the constraint system
    for i = 1:numel(D11)
        A = [A; kron(D21{i}',D12{i})];
        b = [b; -D11{i}(:)];
    end
    
    if (isstruct(structure) && ~isempty(Z) )
        % If Z is not empty it means that the user also imposed a specific structure 
        % on the D matrix of the controller and this structure should be added to the other constraints
        % the constraint is of the type D(i,j) = d_ij
        A = [A; Z];
        b = [b; Z*structure.Dfix(:)];	% here I changed to account for the fixed values of the parameters 
    end
       
    [x,V,w] = solveSystem(A,b);
    
    if (isnan(x))   %no solution
        Dk = nan;
    elseif ~isempty(x)      %unique solution
        Dk = reshape(x,m,p);
    else            %infinity solution, parametrized through V and w
        Dk = [];
    end
    
    if ~isstruct(structure)
        structure = [];
    end
    
    structure.Dk = Dk;
    structure.V = V;
    structure.w = w;
    
    structure.D22 = D22{1};     % save D22 matrix for the final transformation
end

if (isstruct(structure) && isfield(structure,'Dhat') && isempty(structure.Dhat))      
    % If the structure matrix Dhat exists but is empty it means that the
    % user set Dhat = 0
    % If Dk or V matrix are not set it means that we are not in the case of
    % H2 optimization. To use the existing structure we just set Dk=0 as if we had 
    % an unique solution to the constraints system above
    if (~isfield(structure,'Dk') || ~isfield(structure,'V') || ... 
                (isempty(structure.Dk) && isempty(structure.V)) )
        structure.Dk = structure.Dfix; %zeros(m,p);  here I changed to account for the new functionnality
    end
end

%Dk, V and w in the structure include the information from Dhat

structure_out = structure;
