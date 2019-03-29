function [var_bout] = findVarBout(boutIndices, var_temp, dim) 

[n, m] = size(boutIndices);

if dim == 3
    var_bout = nan(n, m, 2);
    for i = 1 : size(boutIndices,1)
        % if 0 in bout indices
        %f = find(ind_b_a_bout(i,:)==0,1)-1;
        % if NaNs
        f = find(isnan(boutIndices(i,:)), 1) - 1;
        if ~isempty(f)
            var_bout(i,1:f,:) = var_temp( boutIndices(i,1:f), :, i);
        else
            var_bout(i,:,:) = var_temp(boutIndices(i,:),:,i);
        end
    end
    
elseif dim == 2
    var_bout = nan(n, m);
    for i = 1:size(boutIndices,1)
        % if 0 in bout indices
        %f = find(ind_b_a_bout(i,:)==0,1)-1;
        % if NaNs
        f = find(isnan(boutIndices(i,:)),1)-1;
        if isempty(f) == 0
            var_bout(i,1:f) = var_temp(i,boutIndices(i,1:f));
        else
            var_bout(i,:) = var_temp(i,boutIndices(i,:));
        end
    end
end