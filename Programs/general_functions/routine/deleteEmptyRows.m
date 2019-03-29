function [Mresized, rowsToDel] = deleteEmptyRows(M, N)

% ::: if 1 input :::
% input : matrix M
%         empty rows = NaNs only
% output : M resized

% ::: if 2 inputs :::
% input : matrix M and N, same size
%         empty rows of N will be deleted in M
% output : M resized
switch nargin
    case 1
        allNaNs = size(M, 2) - sum(isnan(M),2);
        rowsToDel =  find(allNaNs == 0) ;
        M(rowsToDel, :) = [];
        Mresized = M;
        
    case 2
        allNaNs = size(N, 2) - sum(isnan(N),2);
        rowsToDel =  find(allNaNs == 0) ;
        M(rowsToDel, :) = [];
        Mresized = M;
end

end