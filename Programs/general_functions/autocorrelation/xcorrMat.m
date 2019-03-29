function [ACnorm] = xcorrMat (V)

new_start_after_nan = diff(isnan(V),1,2);
[newseq_row, newseq_col] = find(new_start_after_nan == -1, 1);
while ~isempty(newseq_row)
    newseq = V(newseq_row, newseq_col+1:end);
    addnans = size(V,2)-length(newseq);
    V(newseq_row, newseq_col+1:end) = nan(1,length(newseq));
    newseq = [newseq nan(1,addnans)];
    V = [V; newseq];
    new_start_after_nan = diff(isnan(V),1,2);
    [newseq_row, newseq_col] = find(new_start_after_nan == -1, 1);
end

ac = NaN(size(V));
for i = 1 : size(V,1)
    sequence = V(i,:);
    firstnan = find(isnan(sequence),1);
    if firstnan ~= 1
        sequence = sequence(1:firstnan-1);
        xc = xcorr(sequence, 'biased');
        ac(i,1:firstnan-1) = xc( end - length(sequence)+1 : end);
    else
        continue
    end
end

acnorm = ac./ac(:,1);

ACnorm = nanmean(acnorm,1);