function [Astr,Afix,Ano] = getStructureMatrix(Astruct,oldStruct)

%given a structure matrix we construct the matrix that extracts from a matrix M
% the elements to be optimized xM = Astr*M equivalent to M(Astruct) in old version
% Astruct the matrix defining the structure 
% oldStruct is 1 for the or the old format with logicals and 0 for the new format with NaN
% Ast - the new structure matrix obtained by selecting from the I matrix
% all the elements that correspond to parameters that must be optimized
% Afix - matrix of the same size as Astruct containing the values for the
% fixed parameters and 0s for the ones to be optimized
% Ano specifies which elements should not be optimized 

[n,m] = size(Astruct);
npars = n*m;

I = eye(n*m);
vect = Astruct(:);
Astr = []; Ano = []; Afix = [];

if (oldStruct == 1) % old type of structure, Afix remains only zeros
    Afix = zeros(n,m);  % old stucture only allows to fix elements to 0
    for i = 1:npars
        if vect(i) ~= 0
            Astr = [Astr; I(i,:)];
        else
            Ano = [Ano; I(i,:)];        % selects which elements should be fixed
        end
    end
else   % must compare with NaNs
    Afix = vect;
    for i = 1:npars
        if isnan(vect(i))       % parameter that should be optimized
            Astr = [Astr; I(i,:)];
            Afix(i) = 0;
        else
            Ano = [Ano; I(i,:)];
        end
    end
    Afix = reshape(Afix,n,m);
end
