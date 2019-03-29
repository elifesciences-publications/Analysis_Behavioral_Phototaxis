function[] = bias_distribution_individuals(different_fish, dXi, FishID, sequencesperfish, fishdXmean, fishdXstd, save_fig)

meandXperFishi_onseq = NaN(length(different_fish), 1);
stddXperFishi_onseq = NaN(length(different_fish), 1);

binwidth = 0.1;

%***
fig1 = figure;
for i = 1 : length(different_fish)
    rowsOfsequencesOneFish = find(FishID == i);
    dxi = dXi(rowsOfsequencesOneFish,:);
    meandXperSeqi = nanmean(dxi, 2);
    
    meandXperFishi_onseq(i) = nanmean(meandXperSeqi);
    stddXperFishi_onseq(i) = nanstd(meandXperSeqi);
    
    [yhist, edges]  = histcounts(dxi, 30, 'Normalization', 'probability');
    x_hist=edges(1:end-1)+binwidth/2;
    yhist=yhist/(sum(yhist)*binwidth);
%     plot(x_hist,smooth(yhist,7), ...
%         'Color', [i/length(different_fish)/2 i/length(different_fish)/2 i/length(different_fish)/2]) 
%     hold on
    %waitforbuttonpress
    
    % -- individual sequence/fish biases --
    plot(i*ones(length(meandXperSeqi)), meandXperSeqi, 'o')
    hold on
end
axes1 = fig1.Children;
axes1.Title.String = ['probability distribution of \delta\theta for ' num2str(i) ' fish'];
axes1.XLim = [- pi pi];
axes1.XAxis.Label.String = '\delta\theta';
axes1.XAxis.TickValues = -pi:pi/2:pi;
axes1.XAxis.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
axes1.YAxis.Label.String = 'pdf';
axes1.FontSize = 14;
axes1.FontName = 'Times New Roman';

%*** 
fig2 = figure;

subplot(3,1,1)
errorbar(fishdXmean, meandXperFishi_onseq, stddXperFishi_onseq,...
    'o','LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.2 0.2 0.2])
hold on
errorbar(fishdXmean, meandXperFishi_onseq, stddXperFishi_onseq./sequencesperfish,...
    'o','LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.8 0.2 0.2])
plot([-0.2 0.2], [-0.2 0.2], '--', 'Color', [0.7 0 0.2])
xlabel('<\delta\theta> on all bouts of individual fish (rad)')
ylabel('<\delta\theta> on sequences of individual fish')
legend('mean and std', 'mean and sem', 'y=x')

subplot(3,1,2)
[yhist, edges]  = histcounts(fishdXmean, 3*round(sqrt(length(fishdXmean))), 'Normalization', 'probability');
x_hist=edges(1:end-1)+binwidth/2;
x_hist = [x_hist(1)-binwidth/2, x_hist, x_hist(end)+binwidth/2];
yhist=yhist/(sum(yhist)*binwidth);
yhist = [0, yhist, 0];
plot(x_hist,smooth(yhist,3), 'k', 'Linewidth', 1.5)
text(5, 30, ['Mean = ' num2str(mean(fishdXmean)) ' deg'])
title('Distribution of  A = < \delta X > ')
xlabel('A = <\delta X > (rad)')
text(0.2, 0.4, ['Mean = ' num2str(mean(fishdXmean)) ' rad'])
title('Mean bout bias for each fish distribution')
xlim([-0.5 0.5])

subplot(3,1,3)
[yhist, edges]  = histcounts(fishdXmean./fishdXstd, 3*round(sqrt(length(fishdXmean))), 'Normalization', 'probability');
x_hist=edges(1:end-1)+binwidth/2;
x_hist = [x_hist(1)-binwidth/2, x_hist, x_hist(end)+binwidth/2];
yhist=yhist/(sum(yhist)*binwidth);
yhist = [0, yhist, 0];
plot(x_hist,smooth(yhist,3), 'k', 'Linewidth', 1.5)
title('Distribution of  <\delta\theta> / \surd\delta\theta')
xlabel('<\delta\theta>/\surd\sigma(\delta\theta)')
xlim([-1 1])

fig2.Name = 'Bias during spontaneous swimming';
for i = 1 : size(fig2.Children, 1)
    fig2.Children(i).FontSize = 14;
    fig2.Children(i).FontName = 'Times New Roman';
end

if save_fig
    save_fig_and_svg(directory, name)
end