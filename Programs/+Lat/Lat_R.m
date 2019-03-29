
%% after AnalysisLat

red = [0.7 0 0];

%%
% --- distribution of X between 'bornes' ---

bornes = [3 20];
x = Xfilt(:, bornes(1) : bornes(2));

% pdf
Nbins = 18;
angle = ((1:Nbins) - 0.5) * 2*pi/Nbins;
Rpdf = hist(mod(x(:),2*pi),Nbins)/sum(hist(mod(x(:),2*pi),Nbins));

% R
theta_mean = circ_mean(x(:));
Rproj_mean = circ_r(x(:)).*cos(theta_mean+pi/2);

% text on plot
txtAngle = [pi/2, pi/2];
txtRng = [0.04, 0.03];
names = {['<\theta> = ' num2str(round(wrapToPi(theta_mean-pi/2)*100)/100) ' rad'],...
    ['<R> = ' num2str(round(Rproj_mean*100)/100)]};

%***
fig = figure;
handles.tab2_axes = axes;
polarplot([angle, angle(1)]-pi/2,[Rpdf, Rpdf(1)],...
    'Linewidth', 2, 'Color', [0 0 0]);
hold on
polarplot([0 theta_mean]-pi/2, [0 Rproj_mean/10], 'Color', red, 'Linewidth', 1.5)
text(txtRng.*cos(txtAngle),txtRng.*sin(txtAngle), names,...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize', 14, 'FontName', 'Times New Roman')

set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','clockwise')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.ThetaAxisUnits = 'radians';
ax.RAxisLocation = -3*pi/4;    
thetaticks((0 : pi/4 : 2*pi))
thetaticklabels({'0', '', '\pi/2', '', '\pi', '', '-\pi/2'})
rtickangle(0)
title(['bouts # ' num2str(bornes(1)) ' to ' num2str(bornes(2))])

%% save figure
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/FiguresSerious/';
name = ['lateralized_Xdistribution3_11'];
saveas(fig,[figpath name],'fig')
disp('fig saved')

%%
% --- Resultant vector length ---
meanXovertime = circ_mean(Xfilt);
R_meant = circ_r(Xfilt);
Rproj_meant = R_meant.*cos(meanXovertime + pi/2);

%***
fig = figure;
plot(R_meant(1:50), 'Linewidth', 2)
hold on
plot(Rproj_meant(1:50), 'Linewidth', 2)
xlabel('bout #')
legend('R : mean resultant vector length', 'R.cos(<X>)')
xlim([0 20])
grid on
ax = gca;
ax.FontSize = 14;

%% R per fish
meanX_perfish = nan(length(unique(FishN)), 1);
Rperfish = nan(length(unique(FishN)), 1);
for i = unique(FishN)'
	fish = find(FishN == i);
    xfish = XLat(fish, bornes(1):bornes(2) );
    if sum(isnan(xfish(:))) == numel(xfish)
        continue
    end
    meanX_perfish(i) = circ_mean(xfish(:));
    Rperfish(i) = circ_r(xfish(:));
end
Rproj_perfish = Rperfish.*cos(meanX_perfish+pi/2);

Nbins = 9;
rbins_perfish = ((-Nbins/2:Nbins/2-1) + 0.5) * range(Rproj_perfish(:))/Nbins;
Rpdf_perfish = hist(Rproj_perfish,Nbins)./sum(hist(Rproj_perfish,Nbins));

%***
figure;
bar(sort(Rproj_perfish))
sum(Rproj_perfish>0)/length(Rproj_perfish)
text(5, 0.6, ['Rsource>0 ' num2str(sum(Rproj_perfish>0)/length(Rproj_perfish)*100) '% of the fish'])
xlabel('fish #')
ylabel('Rsource = R.cos(<X>-pi)')

%***
fig = figure;
plot(rbins_perfish, Rpdf_perfish,'k', 'Linewidth', 1.5)
hold on
plot([0 0], [0 0.35], 'k')
text(0.6, 0.3, 'N = 38', 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel('pdf')
xlabel('mean R for individual fish')
title('Resultant vector of <X> relative to source, for individual fish')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
%% save figure
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/FiguresSerious/';
name = ['lateralized_RperfishFish'];
saveas(fig,[figpath name],'fig')
disp('fig saved')
%% R per fish (2)

% --- ksdensity ---
[Rpdf_perfish_ks, rbins_perfish_ks] = ksdensity(Rproj_perfish,...
    min(Rproj_perfish):0.1:max(Rproj_perfish), 'Bandwidth', 0.1, 'Kernel', 'normal');

%*** 
fig = figure;
plot(rbins_perfish_ks, Rpdf_perfish_ks,'k', 'Linewidth', 1.5)
hold on
plot([0 0], [0 max(Rpdf_perfish_ks)*1.1], 'k')
text(0.6, max(Rpdf_perfish_ks)*0.8, 'N = 38', 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel('pdf')
title('Resultant vector of <X> relative to source, for individual fish')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%% save figure
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/FiguresSerious/';
name = ['lateralized_RpdfFish_ksdensity'];
saveas(fig,[figpath name],'fig')
disp('fig saved')