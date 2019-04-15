%% MSD

% Mean squared displacement analysis
% on spontaneous swimming data
% .........................................................................

%% Check coefficient corresponding to individual P_switch...
p_flip_perfish = NaN(size(aciPF,1),1);
data_to_fit_on = 2;

for i = 1 : size(aciPF,1)
    int = [0:9];
    acf = aciPF(i,int+1);
    [autocorr_fit]  = AutocorrelationFit(data_to_fit_on, int, acf, pturn, wturn, wfor, 0);
    p_flip_perfish(i) = autocorr_fit.p_flip;
end
root_path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/';
folder = 'SimuData';
load([root_path folder filesep 'SimuPswitch.mat'], 'SimuPswitch')
psfit = SimuPswitch.fit;

correspondingcoeff = psfit(p_flip_perfish);
%% MSD per fish :
% Compare D obtained from fitting MSD
% to D = var*int(AC)

di = NaN(length(different_fish),1);
<<<<<<< HEAD
=======
max = 1;
>>>>>>> master
MSDperfishi = NaN(length(different_fish), size(Xi,2));
for fish = 1 : length(different_fish)
    rowsOfsequencesOneFish = find(FishID == fish);
    Xifish = Xi(rowsOfsequencesOneFish, :);
    fb = fishBias(rowsOfsequencesOneFish);
    MSDperfishi(fish,:) = msdX0shuffled(Xifish, fb); 
    ll = 3;
    ul = 15;
    linfiti = polyfit((ll:ul), MSDperfishi(i,ll:ul),1);
    di(fish) = linfiti(1);
end

meanMSDperfishi = nanmean(MSDperfishi,1);
stdMSDperfishi = nanstd(MSDperfishi,1,1);
nMSDperfishi = length(different_fish) - sum(isnan(MSDperfishi),1);

%%
% --- Diffusion coeff : SAC.std VS MSD
fig = figure;
plot(0.5 * 0.7 * di, 0.5 * fishdXstd.^2.*nansum(aciPF,2), 'sq',...
    'MarkerEdgeColor', 'k', 'Markersize', 8, 'MarkerFaceColor', [0.2 0.2 0.2])
hold on
plot([0.001 1], [0.001 1], '--', 'Color', [0.8 0.8 0.8])
plot(0.5 * correspondingcoeff.*di, 0.5 * fishdXstd.^2.*nansum(aciPF,2), 'sq',...
    'MarkerEdgeColor', 'k', 'Markersize', 8, 'MarkerFaceColor', 'r')

xlabel('D from mean squared displacement (rad^2/bout)')
ylabel('D variance and autocorrelation (rad^2/bout)')
legend([num2str(length(different_fish)) ' fishes'], 'y=x')
ax = gca;
ax.FontSize = 14;
ax.YScale = 'log';
ax.XScale = 'log';

%% MSD on whole block of data

x = Xi;
%x(TimeBouti>100)=NaN;
fb = fishBias;
[MSD_shuff,sem] = msdX0shuffled(x, fb);

% linear fit of MSD
ll = 8;
ul = 20;
linfitMSD = polyfit((ll:ul), MSD_shuff(ll:ul),1);

% --- add simu and analyt. ---
[~, ~, ~, t, ~, MSD] = dx_vs_dx(wturn, wfor, pturn, p_flip);

% ***
bforplot = 0:20;
fig = figure;
hold on
% errorbar(bforplot, MSD_shuff(bforplot+1), sem(bforplot+1),...
%  'Linewidth', 2, 'Color', [0.2 0.2 0.2], 'DisplayName', 'MSD all data')
% plot(bforplot, (bforplot+1)*linfitMSD(1)+linfitMSD(2),...
%  'g--', 'Linewidth', 1.5, 'DisplayName', ['linear fit : a = ' num2str(linfitMSD(1)) ' bouts ' num2str(ll) '-' num2str(ul)]) % linear fit
errorbar(bforplot, meanMSDperfishi(bforplot+1), stdMSDperfishi(bforplot+1)./sqrt(nMSDperfishi(bforplot+1)),...
    'Linewidth', 2, 'Color', colour(2,:), 'DisplayName', '<MSD>_{fish}')
plot(bforplot, bforplot*MSD_shuff(2),...
     '--','Color', colour(4,:), 'DisplayName', 'extension of 1st bout')
 
 linfitMSDpf = polyfit((ll:ul), meanMSDperfishi(ll:ul),1)

hold on
plot(bforplot,MSD(bforplot+1),...
 '.-', 'Linewidth', 1.5, 'DisplayName', 'Analytical solution', 'Color', colour(3,:));
% plot(bforplot,MSDsim(bforplot+1),...
%  '--', 'Linewidth', 1.5, 'DisplayName', 'Simulation');

xlim([0 20])
<<<<<<< HEAD
ylabel('< (\theta_n-\theta_{n_0}) ^2>')
=======
ylabel('< (\delta\theta_n-\delta\theta_{n_0}) ^2>')
>>>>>>> master
xlabel('bout #')
ax =gca;
ax.FontSize= 14;
ax.FontName = 'Times New Roman';

legend


