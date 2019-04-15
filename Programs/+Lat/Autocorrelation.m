%% Extract autocorrelative parameters on bouts from Lateralized (stereo-phototaxis) data
% for lag = 1
% and ACF (lag = 0 : inf)

%% ::: Autocorrelation in reinforcement/inhibition conditions :::
<<<<<<< HEAD
dX = dXLssbiais;
=======
dX = dXl;
>>>>>>> master
C = DIlr;

binsdXdX = 9;
binspflip = 7;
<<<<<<< HEAD
Lat.auto_co_reinforcement(dX, dXissbiais, C, binsdXdX, binspflip, 1, 'absolute')
=======
Lat.a.auto_co_reinforcement(dX, dXissbiais, C, binsdXdX, binspflip, 1)

Lat.a.auto_co_reinforcement_abscontrast(dXl, dXi, C, binspflip+3, 1)
>>>>>>> master

%% Autocorrelation in coherence/conflict situations : NOT INTERESTING
dC = diff(DIlr,1,2);
Cn = DIlr(:, 1:end-2)/max(abs(DIlr(:)));
Cnp1 = DIlr(:, 2:end-1)/max(abs(DIlr(:)));
dXn = dXl(:,1:end-1);
dCn = dC(:,1:end-1);
dXnp1 = dXl(:,2:end);
dCnp1 = dC(:,2:end);
ACn_np1 = dXn.*dXnp1./abs(dXn.*dXnp1);

coherence = find( (Cn > 0 & dXn > 0 & dCn < 0) | (Cn < 0 & dXn < 0 & dCn < 0));
conflict = find( (Cn > 0 & dXn > 0 & dCn > 0) | (Cn < 0 & dXn < 0 & dCn > 0));

% --- dXn+1 VS dXn ---
[binvalsA, elts_per_binA, v2binA] = BinsWithEqualNbofElements(dXn(coherence), dXnp1(coherence), 12, 18);
[binvalsB, elts_per_binB, v2binB] = BinsWithEqualNbofElements(dXn(conflict), dXnp1(conflict), 12, 18);
mV2A = nanmean(v2binA,2);
mV2B = nanmean(v2binB,2);

%*** 
figure;

errorbar(binvalsA, mV2A, std(v2binA,1,2)/sqrt(elts_per_binA), 'Linewidth', 2)
hold on
errorbar(binvalsB, mV2B, std(v2binB,1,2)/sqrt(elts_per_binB), 'Linewidth', 2)
xlabel('dXn')
ylabel('dXn+1')

% --- <dXn+1 dXn> ---

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements((Cnp1), ACn_np1, 6, 9);
[binvalsA, elts_per_binA, v2binA] = BinsWithEqualNbofElements((Cnp1(coherence)), ACn_np1(coherence), 9, 12);
[binvalsB, elts_per_binB, v2binB] = BinsWithEqualNbofElements((Cnp1(conflict)), ACn_np1(conflict), 9, 12);
mV2 = nanmean(v2bin,2);
mV2A = nanmean(v2binA,2);
mV2B = nanmean(v2binB,2);

%***
figure
hold on
errorbar(binvals, mV2, std(v2bin,1,2)/sqrt(elts_per_bin),...
    'k', 'Linewidth', 2, 'DisplayName', 'all')
errorbar(binvalsA, mV2A, std(v2binA,1,2)/sqrt(elts_per_binA),...
    'Linewidth', 2, 'DisplayName', 'coherence')
errorbar(binvalsB, mV2B, std(v2binB,1,2)/sqrt(elts_per_binB),...
    'Linewidth', 2, 'DisplayName', 'conflict')
%plot([binvalsA(1) binvalsA(end)], [0.2793 0.2793],...
%    '--', 'DisplayName', 'baseline')
legend
xlabel('|I_L-I_R|')
ylabel({'\delta\theta_n*\delta\theta_n_+_1', '---', '|(\delta\theta_n*\delta\theta_n_+_1)|'});
