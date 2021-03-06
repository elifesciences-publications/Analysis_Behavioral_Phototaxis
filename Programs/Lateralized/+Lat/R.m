
colour = colour_palette(0, 3);
%% R in time
x = XLat;
% --- Resultant vector length ---
meanXovertime = circ_mean(x);
Rcirc = circ_r(x);
Rproj = Rcirc.*cos(meanXovertime+pi/2);

%***
fig = figure;
plot(Rcirc(1:40),...
    'DisplayName','R : mean resultant vector length', 'Linewidth', 2, 'Color', colour(2,:))
hold on
plot(Rproj(1:40),...
    'DisplayName', 'R.cos(<\theta>_{bout})', 'Linewidth', 2, 'Color', colour(1,:))
xlabel('bout #')
yticks([-0.1:0.1:0.3])
legend
ax = gca;
ax.FontSize = 14;

%% R per fish
meanX_perfish = nan(length(unique(FishID)), 1);
Rperfish = nan(length(unique(FishID)), 1);
bornes = [2 17];
r=[];
f=[];
for i = unique(FishID)'
	fish = find(FishID == i);
    xfish = XLat(fish, bornes(1):bornes(2) );
    if sum(isnan(xfish(:))) == numel(xfish)
        continue
    end
    meanX_perfish(i) = circ_mean(xfish(:));
    Rperfish(i) = circ_r(xfish(:));
    r = [r; cos(circ_mean(xfish')+pi/2)'.*circ_r(xfish')'];
    f = [f; ones(length(fish),1)*i];
%     polarhistogram(xfish(:))
%     i
%     waitforbuttonpress
end
Rproj_perfish = Rperfish.*cos(meanX_perfish+pi/2);

Nbins = 20;
rbins_perfish = ((-Nbins/2:Nbins/2-1) + 0.5) * range(Rproj_perfish(:))/Nbins;
Rpdf_perfish = hist(Rproj_perfish,Nbins)./sum(hist(Rproj_perfish,Nbins));

%***
figure
for i = unique(FishID)'
    plot(smooth(r(f==i)))
    hold on
end
xlabel('sequence #')
ylabel('R per sequence')

%***
fig1 = figure;

barh(sort(Rproj_perfish), 'FaceColor', colour(1,:))
sum(Rproj_perfish>0)/length(Rproj_perfish)
text(5, 0.6, ['Rsource>0 ' num2str(sum(Rproj_perfish>0)/length(Rproj_perfish)*100) '% of the fish'])
ylabel('fish #')
xlabel('Rsource = R.cos(<\theta>)')
xlim([-1 1])
title(['between bouts ' num2str(bornes)])
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%***
fig2 = figure;
plot(rbins_perfish, smooth(Rpdf_perfish),'k', 'Linewidth', 1.5)
hold on
plot([0 0], [0 0.15], 'k')
text(0.6, 0.3, ['N = ' num2str(length(unique(FishID)))],...
    'FontSize', 14, 'FontName', 'Times New Roman')
ylabel('pdf')
xlabel('mean R for individual fish')
xlim([-1 1])
title('Resultant vector of <\theta> relative to source, for individual fish')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%% save figure
if 0
    name = ['lateralized_RperfishFish'];
    saveas(fig,[figpath name],'svg')
    saveas(fig,[figpath name],'fig')
    disp('fig and svg saved')
end
%% R per fish
meanXperfish = NaN(length(unique(FishID)),1);
Rcircperfish = NaN(length(unique(FishID)),1);
nbouts = NaN(length(unique(FishID)),1);
stdRonseqperfish = NaN(length(unique(FishID)),1);
meanRonseqperfish = NaN(length(unique(FishID)),1);
seqperfish = NaN(length(unique(FishID)),1);
h_mean = NaN(length(unique(FishID)),1);
pval_nonunif = NaN(length(unique(FishID)),1);

for i = unique(FishID)'
    fish = find(FishID == i);
    xfish = XLat(fish,2:17);
    meanRonseqperfish(i) = nanmean(circ_r(xfish'));
    stdRonseqperfish(i) = nanstd(circ_r(xfish'));
    seqperfish(i) = length(xfish(:,1)) - sum(isnan(xfish(:,1)));
    xfish = xfish(:);
    xfish(isnan(xfish))=[];
    meanXperfish(i) = circ_mean(xfish);
    Rcircperfish(i) = circ_r(xfish);
    nbouts(i) = numel(xfish);
    
    [h_mean(i), ~] = circ_mtest(xfish, pi/2);
    [pval_nonunif(i), ~] = circ_rtest(xfish);
end
Rprojperfish = Rcircperfish.*cos(meanXperfish+pi/2);
meanRprojonseqperfish = meanRonseqperfish.*cos(meanXperfish+pi/2);

%***
figure
[sorted_angle, sorted_idx] = sort(meanXperfish);
polarwitherrorbar(sorted_angle'+pi/2,Rcircperfish(sorted_idx)',stdRonseqperfish',colour(1,:))
view([90 -90])

%***
figure
polarplot(meanXperfish+pi/2, Rcircperfish, 'ko', 'MarkerFaceColor', colour(1,:), 'MarkerSize', 10)
hold on
polarplot(meanXperfish+pi/2, meanRonseqperfish, 'ko', 'MarkerFaceColor', colour(1,:)*0.1, 'MarkerSize', 10)
set(gca,'ThetaZeroLocation','top',...
    'ThetaDir','clockwise')
ax=gca;
ax.ThetaAxisUnits = 'radians';
ax.RAxisLocation = pi;
ax.RMinorGrid = 'on';
rticks(0:0.5:1)
thetaticks((0 : pi/4 : 2*pi))
thetaticklabels({'0', '', '\pi/2', '', '\pi', '', '-\pi/2'})

ax.FontName = 'Times New Roman';
ax.FontSize = 16;

%***
figure
plot(wrapToPi(meanXperfish+pi/2), Rprojperfish,...
    'ko', 'MarkerFaceColor', colour(1,:), 'MarkerSize', 10, 'DisplayName', 'Mean on all bouts/individual')
hold on
errorbar(wrapToPi(meanXperfish+pi/2),meanRprojonseqperfish, stdRonseqperfish./sqrt(seqperfish),...
    'Color', colour(2,:), 'Marker', 'o', 'MarkerFaceColor', colour(1,:), 'Linewidth', 2, 'LineStyle', 'none',...
    'DisplayName', 'Mean & SEM on different sequences of individuals')
plot([-pi pi], [0 0], '--k', 'HandleVisibility', 'off');
plot([0 0], [-1 1], '--k', 'HandleVisibility', 'off');
xlim([-pi pi])
xticks([-pi,-pi/2, 0,pi/2, pi])
xticklabels({'-\pi','-\pi/2', '0','\pi/2', '\pi'})
xlabel('<\theta>_f_i_s_h')
ylabel('R_f_i_s_h')
ax=gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 18;

% Combine non-uniformity with mean direction
h_nonunif = 1-pval_nonunif;
h_nonunif(h_nonunif>=0.99) = 1;
h_nonunif(h_nonunif<0.99) = 0;
signif_phototaxis = sum(h_mean.*h_nonunif')/length(h_mean);
