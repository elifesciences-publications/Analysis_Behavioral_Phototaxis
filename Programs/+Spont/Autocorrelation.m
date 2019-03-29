% Autocorrelation analysis
% of spontaneous swimming data

%% autocorrelation (1/2) 
[ACinorm, sem] = xcorrMatrixRows (dXissbiais);

% fish mean autocorrelation
[ACinormpf, sempf, acPF] = xcorrMatrixRowsPerFish (dXissbiais, FishID);

%% autocorrelation (2/2) : lag1
dT = 1; % lag
Spont.a.autocorr_lag_dt(dXissbiais,dT);

%% Fit autocorrelation function to get p_flip

data_to_fit_on = 2; %2 AC function; 1 AC lag
[autocorr_fit]  = AutocorrelationFit(data_to_fit_on, [0:20], ACinormpf(1:21), pturn, wturn, wfor, 1);
p_flip = autocorr_fit.p_flip;
autocorr_fit

%% autocorrelation with simu (1/2) : function 

% --- from analytical & simu --- 
[x, C1, convC1err, t_a, Cor, MSD] = dx_vs_dx(wturn, wfor, pturn, p_flip);
%[DX1, DY1, t_sim, ac, MSDsim] = simulate_markov(wturn, wfor, pturn, p_flip);

% inputs 4 & 5 : analytical solution (optional)
%        6 & 7 : simulation (optional)
fig = Spont.a.autocorr_function_plots(dXissbiais, FishID, t_a, Cor); %, t_sim, ac

txt = text(5, p_flip, {['p_f_l_i_p = ' num2str(p_flip, 3)]}) ;
txt.FontSize = 14;
txt.FontName = 'Times New Roman';

save_fig=0;
if save_fig
    save_fig_and_svg(fig, 'final/spont', 'ACF');
    save_fig_png(fig, 'final/spont', 'ACF');
end
%% autocorrelation with simu (2/2) : lag1
[colour] = colour_palette(0,1);
dT = 1; % lag
fig = Spont.a.autocorr_lag_dt(dXissbiais,dT);

[errfit] = theta_measure_error_estimation(x);
%convDY1err = conv(DY1,errfit, 'same')/sum(errfit);

% --- add analyt. and simu ---
hold on
 plot(x,convC1err,...
     'DisplayName', 'C(1)\otimes err', 'Linewidth', 1.5, 'Color', colour(4,:));
%plot(x,smooth(C1),...
%    'DisplayName', 'C(1) : analytical solution', 'Linewidth', 1, 'Color', colour(1,:));
% plot(DX1,DY1, '--', ...
%     'DisplayName', 'simulation', 'Linewidth', 1.5, 'Color', colour_simu(1,:));
%plot(DX1,convDY1err, '-.', 'Linewidth', 1.5, 'DisplayName', 'simu convolved');

legend
%%
% --- per fish ---
[~, ~, aciPF] = xcorrMatrixRowsPerFish (dXissbiais, FishID);

% ***
figure;
title('check AC per fish')
plot(aciPF')

%% Calculations of pflip from c(1)

dXn = dXissbiais(:,1:end-1);
dXnp1 = dXissbiais(:,2:end);
pdXndXnp1 = dXn.*dXnp1;
mdXndXnp1 = nanmean(pdXndXnp1(:));
pabsdXndXnp1 = abs(dXn).*abs(dXnp1);
mabsdXndXnp1 = nanmean(pabsdXndXnp1(:));

sem = (nanstd(pdXndXnp1(:))+1/nanstd(pabsdXndXnp1(:)))/sqrt(numel(dXn(:))-sum(isnan(dXn(:))))

pflip1 = 0.5*(1 - mdXndXnp1/(2/pi*pturn^2*wturn^2))

corrfact = (1 - 2*(1-pturn)*wfor/(pturn*wturn));
%corrfact = 0.72
pflip2 = 0.5*(1 - mdXndXnp1/(corrfact*mabsdXndXnp1))

cf = mdXndXnp1/((1-2*p_flip)*mabsdXndXnp1)

%% with inter-bout interval
dXn = dXissbiais(:,1:end-1);
dXnp1 = dXissbiais(:,2:end);
pdXndXnp1 = dXn.*dXnp1;
pabsdXndXnp1 = abs(dXn).*abs(dXnp1);

var1 = IBIi(:,1:end-1);
var2 = pdXndXnp1./pabsdXndXnp1;
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, 12, 20);
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
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, 20, 30);
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