function[XLat,Xfilt,Xlab, FishN, Advancement, TimeBout, xCoord, yCoord, sequencesperfish, boutsperfish]...
    = delete_f_not_enough_bouts(El, threshold)
%%
% --- get variables from structure ---
XLat = El.AngleSource;
Xlab = El.AngleLab;
Xfilt = El.AngleSourceFiltered;
FishN = El.FishN;
Advancement = El.R;
TimeBout = El.TimeBout;
xCoord = El.xCoord;
yCoord = El.yCoord;

%%
% fish identity
f = abs(diff(FishN));
f(f>0) = 1;
initial_FishID = 1+[0 ; cumsum(f)];
different_fish = unique(initial_FishID);

% --- delete fishes with only one sequence or too few bouts ---
boutsperfish = NaN(length(different_fish), 1);
sequencesperfish = NaN(length(different_fish), 1);
for i = 1 : different_fish(end)
    fishseqs = find(initial_FishID == i);
    xi = XLat(fishseqs,:);
    boutsperfish(i) = numel(xi) - sum(isnan(xi(:)));
    sequencesperfish(i) = length(fishseqs);
end

fish_w_too_few_bouts = find(boutsperfish < threshold);
seqs_to_del = [];
for i = 1 : length(fish_w_too_few_bouts)
    seqs_to_del = [seqs_to_del ; find(initial_FishID == fish_w_too_few_bouts(i))];
end
    
XLat(seqs_to_del, :) = [];
FishN(seqs_to_del) = [];
Advancement(seqs_to_del, :) = [];
TimeBout(seqs_to_del, :) = [];
xCoord(seqs_to_del, :) = [];
yCoord(seqs_to_del, :) = [];

sequencesperfish(fish_w_too_few_bouts) = [];
boutsperfish(fish_w_too_few_bouts) = [];

% NEW fish identity
f = abs(diff(FishN));
f(f>0) = 1;
FishN = 1+[0 ; cumsum(f)];

%% clean up

% delete unnecessary nan columns %
nancolumns = find(size(XLat,1) == sum(isnan(XLat),1));
consec = find(diff(nancolumns)==1);

XLat(:,nancolumns(consec))=[];
Xfilt(:,nancolumns(consec))=[];
Xlab(:,nancolumns(consec))=[];
Advancement(:,nancolumns(consec))=[];
TimeBout(:,nancolumns(consec))=[];
xCoord(:,nancolumns(consec))=[];
yCoord(:,nancolumns(consec))=[];

