function[Xi, FishNi, Ri, TimeBouti, xCoordi, yCoordi, FishID, sequencesperfish, boutsperfish] = delete_f_not_enough_bouts(Es, threshold)

%%
% --- get variables from structure ---
Xi = Es.Angle;
FishNi = Es.FishN;
Ri = Es.R;
TimeBouti = Es.TimeBout;
xCoordi = Es.Coordinates.x;
yCoordi = Es.Coordinates.y;

%%
% fish identity
f = abs(diff(FishNi));
f(f>0) = 1;
initial_FishID = 1+[0 ; cumsum(f)];
different_fish = unique(initial_FishID);

% --- delete fishes with only one sequence or too few bouts ---
boutsperfish = NaN(length(different_fish), 1);
sequencesperfish = NaN(length(different_fish), 1);
for i = 1 : different_fish(end)
    fishseqs = find(initial_FishID == i);
    xi = Xi(fishseqs,:);
    boutsperfish(i) = numel(xi) - sum(isnan(xi(:)));
    sequencesperfish(i) = length(fishseqs);
end

fish_w_too_few_bouts = find(boutsperfish < threshold);
seqs_to_del = [];
for i = 1 : length(fish_w_too_few_bouts)
    seqs_to_del = [seqs_to_del ; find(initial_FishID == fish_w_too_few_bouts(i))];
end
    
Xi(seqs_to_del, :) = [];
FishNi(seqs_to_del) = [];
Ri(seqs_to_del, :) = [];
TimeBouti(seqs_to_del, :) = [];
xCoordi(seqs_to_del, :) = [];
yCoordi(seqs_to_del, :) = [];

sequencesperfish(fish_w_too_few_bouts) = [];
boutsperfish(fish_w_too_few_bouts) = [];

% NEW fish identity
f = abs(diff(FishNi));
f(f>0) = 1;
FishID = 1+[0 ; cumsum(f)];

%% clean up

% delete unnecessary nan columns %
nancolumns = find(size(Xi,1) == sum(isnan(Xi),1));
consec = find(diff(nancolumns)==1);

Xi(:,nancolumns(consec))=[];
Ri(:,nancolumns(consec))=[];
TimeBouti(:,nancolumns(consec))=[];
xCoordi(:,nancolumns(consec))=[];
yCoordi(:,nancolumns(consec))=[];

