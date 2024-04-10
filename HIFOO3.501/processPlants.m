function [plantinfo, isstrongstab, ssinput,m,p] = processPlants(plant,ubounds, obj,options)
%Input: plants-- Cell array of plants with possibility of 'K' for strong stabilization
%                in proper internal representation for HIFOO
%       ubounds -- upperbound constraints.  NOTE:  This argument is no longer required
%                  but maintained for backwards compatibility.
%       options-- hifoo options
%               has fields options.optimize and options.constraints 
%               giving indices of plants.
%Output:
%       plantinfo -- cell array of plants, with proper internal format
%                    Each plant, P, has additional fields 
%                    -- P.obj giving objective function
%                    -- (if constraint) P.constraint gives upperbound
%       isstrongstab -- 1 is strong stabilization is selection, 0 otherwise
%       Will throw error if the plants are not of compatible dimensions

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

isstrongstab = 0;
plantinfo = cell(numel(plant),1);
prtlevel = options.prtlevel;
mvals = [];
pvals = [];

ssinput = 0;
for i = 1:numel(plant)
    thisplant = plant{i};
    
    if ischar(thisplant) % plant is a string giving COMPleib name or 'K'
    
    % first determine if it is the controller
     if (strcmp(thisplant,'K') || strcmp(thisplant,'k'))
         isstrongstab = 1;
         P.controller = 1;
         P.A = 0;
         P.B1 = 0;
         P.B2 = 0;
         P.C1 = 0;
         P.C2 = 0;
         P.D11 = 0;
         P.D12 = 0;
         P.D21 = 0;
         P.D22 = 0;
         
         
    % convert from COMPleib format to ours
    else 
        if exist('COMPleib')~= 2
            error('hifoo: when input "plant" is a string "COMPleib" must be installed; see www.compleib.de')
        end 
        [A,B1,B2,C1,C2,D11,D12,D21]=COMPleib(thisplant);
        if isempty(A)
            error('hifoo: input "plant" is a string but is not recognized by COMPleib')
        end
        %pars = [];
            P.A = A;
            P.B1 = B1;
            P.B2 = B2;
            P.C1 = C1;
            P.C2 = C2;
            P.D11 = D11;
            P.D12 = D12;
            P.D21 = D21;
            P.D22 = 0;
        m = size(P.B2,2);
        p = size(P.C2,1);
        mvals = [mvals m];
        pvals = [pvals p];
        P.controller = 0;
       %All compleip examples have D22 = 0;
   end

    elseif strcmp(class(thisplant), 'ss') % SS object
        [A, B, C, D] = ssdata(thisplant);  % standard toolbox conversion
        plantstruct = get(thisplant); % convert to struct
        plantstructIn = plantstruct.InputGroup;
        if ~isfield(plantstructIn, 'U1')|~isfield(plantstructIn, 'U2')
            if prtlevel > 0
                fprintf('hifoo: "plant" is an SS object but cannot find U1 and U2 fields in "InputGroup"\n')
            end
            if (obj(i) == 'h' | obj(i) == 't')
                error('hifoo: cannot optimize H-infinity or H2 performance when input or output performance channels are unspecified') 
            else
                U1 = [];
                U2 = 1:size(B,2);
                if prtlevel > 0
                    fprintf('hifoo: partitioning B = [B1 B2] with B1 empty\n')
                end
            end
        else
            U1 = plantstructIn.U1;
            U2 = plantstructIn.U2;
            if any(sort([U1 U2]) ~= (1:size(B,2)))
                error('hifoo: input "plant" is an SS object but U1 and U2 fields in "InputGroup" do not provide a valid input partitioning')
            else
                if prtlevel > 0
                    fprintf('hifoo: "plant" is an SS object, using input partitioning in "InputGroup"\n')
                end
            end
        end
        plantstructOut = plantstruct.OutputGroup;
        if ~isfield(plantstructOut, 'Y1')|~isfield(plantstructOut, 'Y2')
            if prtlevel > 0
                fprintf('hifoo: "plant" is an SS object but cannot find Y1 and Y2 fields in "OutputGroup"\n')
            end
            if obj(i) == 'h' 
                error('hifoo: cannot optimize H-infinity performance when output performance channel is unspecified')
            else    
                Y1 = [];
                Y2 = 1:size(C,1);
                if prtlevel > 0
                    fprintf('hifoo: partitioning C = [C1; C2] with C1 empty\n')
                end
            end
        else
            Y1 = plantstructOut.Y1;
            Y2 = plantstructOut.Y2;
            if any(sort([Y1 Y2]) ~= (1:size(C,1)))
                error('hifoo: input "plant" is an SS object but Y1 and Y2 fields in "OutputGroup" do not provide a valid output partitioning')
            else
                if prtlevel > 0
                    fprintf('hifoo: "plant" is an SS object, using output partitioning in "OutputGroup"\n')
                end
            end
        end

        P = [];
        P.A = A;
        P.B1 = B(:,U1);
        P.B2 = B(:,U2);
        P.C1 = C(Y1,:);
        P.C2 = C(Y2,:);
        P.D11 = D(Y1,U1);
        P.D12 = D(Y1,U2);
        P.D21 = D(Y2,U1);
        P.D22 = D(Y2,U2);
        n = length(P.A);
        m = size(P.B2,2);
        p = size(P.C2,1);
        mvals = [mvals m];
        pvals = [pvals p];
        ssinput = 1;  % used at end in determining output format
        P.controller = 0;

    elseif isstruct(thisplant) % A, B2 (or B) and C2 (or C) are always required
        if isfield(thisplant,'A')
            P.A = thisplant.A;
        elseif isfield(thisplant, 'a')
            P.A = thisplant.a;
        else
            error('hifoo: input "plant" is a structure but it does not have a field "A"')
        end
        if isempty(P.A)
            error('hifoo: field "A" of input "plant" is empty')
        end
        % field B is a synonym for B2, and C for C2
        if isfield(thisplant,'B2')
            P.B2 = thisplant.B2;
        elseif isfield(thisplant,'B')
            P.B2 = thisplant.B;
        elseif isfield(thisplant,'b2')
            P.B2 = thisplant.b2;
        elseif isfield(thisplant,'b')
            P.B2 = thisplant.b;
        else
            error('hifoo: input "plant" is a structure but it does not have a field "B2" or "B"')
        end
        if isempty(P.B2)
            error('hifoo: field "B2" (or "B") of input "plant" is empty')
        end
        if isfield(thisplant,'C2')
            P.C2 = thisplant.C2;
        elseif isfield(thisplant,'C')
            P.C2 = thisplant.C;
        elseif isfield(plant,'c2')
            P.C2 = thisplant.c2;
        elseif isfield(thisplant,'c')
            P.C2 = thisplant.c;
        else
            error('hifoo: input "plant" is a structure but it does not have a field "C2" or "C"')
        end
        if isempty(P.C2)
            error('hifoo: field "C2" (or "C") of input "plant" is empty')
        end
        
        
        % Convert P.D to proper format.  This field is optional, so no error is thrown.
        if isfield(thisplant,'D')
            P.D22 = thisplant.D;
        elseif isfield(thisplant, 'd')
            P.D22 = thisplant.d
        else
            P.D22 = 0;
        end
        P.controller = 0;
        n = length(P.A);
        m = size(P.B2,2);
        p = size(P.C2,1);
        mvals = [mvals m];
        pvals = [pvals p];
        % if the objective is not 'h', any B1, C1 or Dij provided 
        % should be ignored, so set them to empty
        if (obj(i) ~= 'h' && obj(i) ~= 't') 
            if (isfield(thisplant, 'B1') | isfield(thisplant, 'C1')) & prtlevel > 0
                fprintf('hifoo: objective is not "h" or "t" so ignoring plant.B1, plant.C1, etc\n')
            end
            P.B1 = zeros(n,0);
            P.C1 = zeros(0,n);
            P.D11 = zeros(0,0);
            P.D12 = zeros(0,m);
            P.D21 = zeros(p,0);
        else % process the other fields (only the case of objective 'h')
            if isfield(thisplant,'b1') & ~isfield(thisplant,'B1')
                thisplant.B1 = thisplant.b1;
            end
            if isfield(plant,'c1') & ~isfield(thisplant,'C1')
                thisplant.C1 = thisplant.c1;
            end
            if isfield(thisplant,'d11') & ~isfield(thisplant,'D11')
                thisplant.D11 = thisplant.d11;
            end
            if isfield(thisplant,'d12') & ~isfield(thisplant,'D12')
                thisplant.D12 = thisplant.d12;
            end
            if isfield(thisplant,'d21') & ~isfield(thisplant,'D21')
                thisplant.D21 = thisplant.d21;
            end
            if ~isfield(thisplant,'B1')|~isfield(thisplant,'C1')|... 
             ~isfield(thisplant,'D11')|~isfield(thisplant,'D12')|~isfield(thisplant,'D21')
                if prtlevel > 0 % provide a little additional information
                    fprintf('hifoo: input "plant" is a structure but one of B1, C1, D11, D12 or D21 is missing\n')
                end
                error('hifoo: cannot optimize H-infinity or H2 performance when input or output performance channels are unspecified')
            elseif isempty(thisplant.B1)|isempty(thisplant.C1)|isempty(thisplant.D11)|...
             isempty(thisplant.D12)|isempty(thisplant.D21)
                if prtlevel > 0 % provide a little additional information
                    fprintf('hifoo: input "plant" is a structure but one of B1, C1, D11, D12 or D21 is empty\n')
                end
                error('hifoo: cannot optimize H-infinity or H2 performance when input or output performance channel is missing')
            else
                P.B1 = thisplant.B1;
                P.C1 = thisplant.C1;
                P.D11 = thisplant.D11;
                P.D12 = thisplant.D12;
                P.D21 = thisplant.D21;
                if isfield(thisplant,'D22') || isfield(thisplant,'D') %If NO D22, set it to 0
                    %fprintf('hifoo: D22 != 0.  Using non-linear code');
                    P.D22 = thisplant.D22;
                else
                    P.D22 = 0;
                end
            end
        end
        if diff(size(P.A)) ~= 0
            error('hifoo: field "A" of input "plant" is not square')
        elseif norm(diff([size(P.A,1) size(P.B1,1) size(P.B2,1)])) > 0
            error('hifoo: row dimension mismatch in fields A, B1, B2 of input "plant"')
        elseif norm(diff([size(P.C1,1) size(P.D11,1) size(P.D12,1)])) > 0
            error('hifoo: row dimension mismatch in fields C1, D11, D12 of input "plant"')
        elseif size(P.C2,1) ~= size(P.D21,1)
            error('hifoo: row dimension mismatch in fields C2, D21 of input "plant"')
        elseif norm(diff([size(P.A,2) size(P.C1,2) size(P.C2,2)])) > 0
            error('hifoo: column dimension mismatch in fields A, C1, C2 of input "plant"')
        elseif norm(diff([size(P.B1,2) size(P.D11,2) size(P.D21,2)])) > 0
            error('hifoo: column dimension mismatch in fields B1, D11, D21 of input "plant"')
        elseif size(P.B2,2) ~= size(P.D12,2)
            error('hifoo: column dimension mismatch in fields B2, D12 of input "plant"')
        end
    else
        error('hifoo: input "plant" must be a structure or an SS object or a string')
    end
    
    if isempty(find(options.optimize == i)) 
        constraints = setdiff(1:numel(plantinfo),options.optimize);
        P.constraint = options.constraints(find(constraints == i));
    end
    P.objective = obj(i);
    plantinfo{i} = P;
    clear P;
end


%Verify that dimensions are compatible

if (max(abs(mvals - circshift(mvals,[0,1])))>0) || (max(abs(pvals - circshift(pvals,[0,1]))>0))
    error('hifoo:  Error:  Multiple plants must have compatible dimensions')
end




