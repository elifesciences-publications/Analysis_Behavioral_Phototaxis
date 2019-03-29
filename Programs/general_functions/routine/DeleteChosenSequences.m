function [seqToDel] = DeleteChosenSequences(date, fish)

% sequences to delete

% lateralization experiments

%%
DeleteSequences = ...
    {'18-01-12', 'fish3', 1 ;
    '18-01-16', 'fish1', 4 ;
    '18-01-17', 'fish1', 7 ;
    '18-01-19', 'fish2', 2 ;
    '18-01-23', 'fish2', 7:9 ;
    '18-01-30', 'fish3', 20 ;
    '18-05-16', 'fish1', 1 ;
    '18-05-23', 'fish2', 3 ;
    };

%%
findDelDate = strfind( DeleteSequences(:, 1) , date);
checkDel = 1-cellfun(@isempty, findDelDate);

seqToDel = [];
if sum(checkDel) > 0
    ind = find(checkDel);
    if strcmp(DeleteSequences(ind, 2), fish)
        seqToDel = cell2mat(DeleteSequences(ind,3));
    end
end

end
