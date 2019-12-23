
different_fish = [1:Nexp];
fishBias = 0;%mean(diff(Xi,1,2),2);

MSDperfishi = NaN(length(different_fish), size(thetaSimu,2));
for fish = 1 : length(different_fish)
    Xifish = thetaSimu(fish, :);
    MSDperfishi(fish,:) = msdX0shuffled(Xifish, 0); 
end

meanMSDperfishi = nanmean(MSDperfishi,1);
stdMSDperfishi = nanstd(MSDperfishi,1,1);
nMSDperfishi = length(different_fish) - sum(isnan(MSDperfishi),1);

bforplot = 0:30;
fig = figure;
hold on
% errorbar(bforplot, MSD_shuff(bforplot+1), sem(bforplot+1),...
%  'Linewidth', 2, 'Color', [0.2 0.2 0.2], 'DisplayName', 'MSD all data')
% plot(bforplot, (bforplot+1)*linfitMSD(1)+linfitMSD(2),...
%  'g--', 'Linewidth', 1.5, 'DisplayName', ['linear fit : a = ' num2str(linfitMSD(1)) ' bouts ' num2str(ll) '-' num2str(ul)]) % linear fit
errorbar(bforplot, meanMSDperfishi(bforplot+1), stdMSDperfishi(bforplot+1)./sqrt(nMSDperfishi(bforplot+1)),...
    'Linewidth', 2, 'Color', [0.2 0.4 0.9], 'DisplayName', 'amplitude autocorrelation') %[0.2 0.4 0.9]
ylabel('MSR M_q')
xlabel('number of bouts M_q')
xticks([0:10:30])
ax1 = gca;
ax1.FontSize = 16;
ax1.FontName = 'Times New Roman';

ll = 3;
ul = 15;
linfitMSDpf = polyfit((ll:ul), meanMSDperfishi(ll:ul),1)


%%
dthetaSimu = diff(thetaSimu,1,2);
[a,b, std_ac, elts] = Spont.autocorr_diff_lag_dt_square(dthetaSimu,1);
figure
hold on
shadedErrorBar(b, a, std_ac/sqrt(elts),...
    'lineprops',{'.','Color', [0.2 0.4 0.9]}) %[0.5 0.5 0.5][0.2 0.4 0.9]
ylabel('<\delta\alpha^2_{n+1}>')
ylim([0.5 1.5])
ax2 = gca;
ax2.XScale = 'log'
ax2.FontSize = 16;
ax2.FontName = 'Times New Roman';

legend('no amplitude autocorrelation', 'amplitude autocorrelation')
