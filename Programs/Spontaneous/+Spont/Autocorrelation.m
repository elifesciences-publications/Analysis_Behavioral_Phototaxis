% Autocorrelation analysis
% of spontaneous swimming data

%% autocorrelation (1/2) 
[ACinorm, sem] = xcorrMatrixRows (dXissbiais);

% fish mean autocorrelation
[ACinormpf, sempf, acPF] = xcorrMatrixRowsPerFish (dXissbiais, FishID);

%% autocorrelation (2/2) : lag1 or +
dT = 1; % lag
[fig, aclag] = Spont.autocorr_lag_dt(dXissbiais,dT);

%% Fit autocorrelation function to get p_flip

data_to_fit_on = 2; %2 AC function; 1 AC lag
[autocorr_fit]  = AutocorrelationFit(data_to_fit_on, [0:20], ACinormpf(1:21), pturn, wturn, wfor, 1);
p_flip = autocorr_fit.p_flip;
autocorr_fit

pflip_pf = NaN(size(acPF,1),1);
for i = 1 : size(acPF,1)
    [autocorr_fit]  = AutocorrelationFit(data_to_fit_on, [0:20], acPF(i,1:21), pturn, wturn, wfor, 1);
    
    pflip_pf(i) = autocorr_fit.p_flip;
end

%% autocorrelation with simu (1/2) : function 

% --- from analytical & simu --- 
[x, C1, convC1err, t_a, Cor, MSD] = dx_vs_dx(wturn, wfor, pturn, p_flip);
%[DX1, DY1, t_sim, ac, MSDsim] = simulate_markov(wturn, wfor, pturn, p_flip);

% inputs 4 & 5 : analytical solution (optional)
%        6 & 7 : simulation (optional)
fig = Spont.autocorr_function_plots(dXissbiais, FishID, 'bp', t_a, Cor); %, t_sim, ac

txt = text(5, p_flip, {['p_f_l_i_p = ' num2str(p_flip, 3)]}) ;
txt.FontSize = 14;
txt.FontName = 'Times New Roman';

%% autocorrelation with simu (2/2) : lag1
[colour] = colour_palette(0,3);
dT = 1; % lag
fig = Spont.autocorr_lag_dt(dXissbiais,dT);

[errfit] = theta_measure_error_estimation(x);
%convDY1err = conv(DY1,errfit, 'same')/sum(errfit);

% --- add analyt. and simu ---
hold on
 plot(x,convC1err,...
     'DisplayName', 'C(1)\otimes err', 'Linewidth', 1.5, 'Color', colour(3,:));
%plot(x,smooth(C1),...
%    'DisplayName', 'C(1) : analytical solution', 'Linewidth', 1, 'Color', colour(1,:));
% plot(DX1,DY1, '--', ...
%     'DisplayName', 'simulation', 'Linewidth', 1.5, 'Color', colour_simu(1,:));
%plot(DX1,convDY1err, '-.', 'Linewidth', 1.5, 'DisplayName', 'simu convolved');

xlim([-pi/2-0.1 pi/2+0.1])
xticks([-pi/2, 0, pi/2])
xticklabels({'-\pi/2','0', '\pi/2'})

legend
ax = gca;
ax.Children(2).Color = 'k'
%%
% --- per fish ---
[~, ~, aciPF] = xcorrMatrixRowsPerFish (dXissbiais, FishID);

% ***
figure;
title('check AC per fish')
plot(aciPF')

%% Calculations of pflip from c(1)
% cf ELN - Behavioral analysis - p 34

dXn = dXissbiais(:,1:end-1);
dXnp1 = dXissbiais(:,2:end);
pdXndXnp1 = dXn.*dXnp1;
mdXndXnp1 = nanmean(pdXndXnp1(:));
pabsdXndXnp1 = abs(dXn).*abs(dXnp1);
mabsdXndXnp1 = nanmean(pabsdXndXnp1(:));

sem = (nanstd(pdXndXnp1(:))+1/nanstd(pabsdXndXnp1(:)))/sqrt(numel(dXn(:))-sum(isnan(dXn(:))))

% - 1 -
pflip1 = 1/2*(1 - mdXndXnp1/(2/pi*pturn^2*wturn^2))

% - 2 -
corrfact = 1 - 2*(1-pturn)*wfor/(pturn*wturn);
pflip2 = 0.5*(1 - mdXndXnp1/mabsdXndXnp1/corrfact)

%% Calculations of pflip - trinarized
% cf ELN - Behavioral analysis - p 34

thresh = 0.22;
dXtrin = dXissbiais;
dXtrin(abs(dXtrin)<thresh) = 0;
dXtrin(dXtrin>0) = 1;
dXtrin(dXtrin<0) = -1;

pturn_estimate = nanmean(abs(dXtrin(:)))

dXn = dXtrin(:,1:end-1);
dXnp1 = dXtrin(:,2:end);
pdXndXnp1 = dXn.*dXnp1;
mdXndXnp1 = nanmean(pdXndXnp1(:));
pabsdXndXnp1 = abs(dXn).*abs(dXnp1);
mabsdXndXnp1 = nanmean(pabsdXndXnp1(:));

% - 1 -
pflip3 = 1/2*(1 - mdXndXnp1/(pturn^2))

pflip4 = 1/2*(1 - mdXndXnp1/(pturn_estimate^2))

plot(thresh, pflip4, '.')



%% with inter-bout interval
dXn = dXissbiais(:,1:end-1);
dXnp1 = dXissbiais(:,2:end);
pdXndXnp1 = dXn.*dXnp1;
pabsdXndXnp1 = abs(dXn).*abs(dXnp1);

var1 = IBIi(:,1:end-1);
var2 = pdXndXnp1./pabsdXndXnp1;
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, 9, 12);
mv2 = mean(v2bin,2);
stdv2 = std(v2bin,1,2);

%*** 
figure;
errorbar(binvals, mv2, stdv2/sqrt(elts_per_bin),...
    'Color', colour(2,:), 'LineWidth', 1.5, 'DisplayName', 'all bouts')

% --- only turns ---
dXn = dXissbiais(:,1:end-1);
dXnp1 = dXissbiais(:,2:end);
dXn(abs(dXn)<0.2)=NaN;
dXnp1(abs(dXnp1)<0.2)=NaN;

pdXndXnp1 = dXn.*dXnp1;
pabsdXndXnp1 = abs(dXn).*abs(dXnp1);

var1 = IBIi(:,1:end-1);
var2 = pdXndXnp1./pabsdXndXnp1;
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, 12, 18);
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