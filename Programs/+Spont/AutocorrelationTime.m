
% convert sequences of angles per bouts
% to angles per frames

dangles_per_bouts = dXissbiais;
corresponding_times = TimeBouti(:,1:end-1);

minimal_sampling_rate = 0.1; %seconds

max_t = max(TimeBouti(:));
longest_T = 0 : minimal_sampling_rate : max_t;
length(longest_T)

dAA = [];
TT = [];
for s = 1 : size(corresponding_times,1)
    corresponding_times_seq = corresponding_times(s,:);
    dangles_per_bouts_seq = dangles_per_bouts(s,:);
    dA = [];
    T = [];
    for t = 1 : length(corresponding_times_seq)-1
        ti = corresponding_times_seq(t);
        tf = corresponding_times_seq(t+1);
        time_sequence = ti : minimal_sampling_rate : tf;
        dt = length(time_sequence);
        dangle_per_time_sequence = ones(1,dt)*dangles_per_bouts_seq(t);
        dA = [dA dangle_per_time_sequence];
        T = [T time_sequence];
    end
    length_nans_for_completion = length(longest_T) - length(dA);
    dA = [dA NaN(1,length_nans_for_completion)];
    T = [T NaN(1,length_nans_for_completion)];
    dAA = [dAA; dA];
    TT = [TT; T];
end

plot(TT',dAA')

%%
[ACinorm, sem] = xcorrMatrixRows (dAA);

lags = [1:200];
plot(lags*minimal_sampling_rate, ACinorm(lags))

%% with inter-bout interval
dAA = fillmissing(dAA, 'previous', 2);
TT = fillmissing(TT, 'previous', 2);
dAA(badseq,:) = [];
TT(badseq,:) = [];

dTT = diff(TT, 1, 2);

dXt = dAA(:,1:end-1);
dXtp1 = dAA(:,2:end);
pdXndXnp1 = dXt.*dXtp1;
pabsdXndXnp1 = abs(dXt).*abs(dXtp1);

var1 = dTT;
var2 = pdXndXnp1./pabsdXndXnp1;
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, 40, 50);
mv2 = mean(v2bin,2);
stdv2 = std(v2bin,1,2);

%*** 
figure;
errorbar(binvals, mv2, stdv2/sqrt(elts_per_bin),...
    'Color', colour(2,:), 'LineWidth', 1.5, 'DisplayName', 'all bouts')

% --- only turns ---
dXt = dXissbiais(:,1:end-1);
dXtp1 = dXissbiais(:,2:end);
dXt(abs(dXt)<0.2)=NaN;
dXtp1(abs(dXtp1)<0.2)=NaN;

pdXndXnp1 = dXt.*dXtp1;
pabsdXndXnp1 = abs(dXt).*abs(dXtp1);

var1 = IBIi(:,1:end-1);
var2 = pdXndXnp1./pabsdXndXnp1;
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, 7, 12);
mv2 = mean(v2bin,2);
stdv2 = std(v2bin,1,2);

%*** 
%figure;
hold on
errorbar(binvals, mv2, stdv2/sqrt(elts_per_bin),...
    'Color', colour(3,:), 'LineWidth', 1.5, 'DisplayName', 'turns only |d\theta|>0.2rad')
xlabel('inter-bout interval (sec)')
ylabel('\delta\theta_n\delta\theta_{n+1}/|\delta\theta_n||\delta\theta_{n+1}|')

legend
ax=gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 14;