%% Orientational bias in temporal phototaxis
%  d(theta) = f(theta)
% and 
%  d(theta) = d(illumination)

colour = colour_palette(0,4);
Vart1 = wrapTo2Pi(X(:,1:20))-pi;
Vart2 = dX(:, 1:20); % dX corrected from individual bias

Vart1 = Vart1(:); Vart1(isnan(Vart2))=[];
Vart2 = Vart2(:); Vart2(isnan(Vart2))=[];
for b = 7:3:30;
yl = '<d\theta_n>';
xl = '\theta_n';

[binvals, elts_per_bin, dXX] = BinsWithEqualNbofElements(Vart1, Vart2, b, b+3);
[~, ~, sigma_error] = theta_measure_error_estimation();

%***
%f = figure;
area([-pi pi], [sigma_error/2, sigma_error/2], -sigma_error/2,...
     'FaceColor', colour(4,:), 'FaceAlpha', 0.3, 'LineStyle', 'none', 'ShowBaseLine', 'off',...
     'DisplayName', '\sigma_{\theta error}')
 hold on
errorbar(binvals,  mean(dXX, 2), std(dXX,1,2)/sqrt(elts_per_bin-1),...
    'Linewidth', 1.5, 'Color', colour(3,:),...
    'DisplayName', 'mean bias per bins in all temporal experiments')
 xlim([-pi pi])
 xticks([-pi : pi/2 :pi])
 xticklabels({'-\pi','-\pi/2', '0','\pi/2', '\pi'})
ylim([-0.1 0.1])
end
xlabel(xl)
ylabel(yl)

ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

ylabel('<\delta\theta>')
xlabel('\theta');
%% <dX> VS X per fish
Vart1 = wrapTo2Pi(X(:, 1:end-1));
Vart2 = dX_corr(:, 1:end); % dX corrected from individual bias
Vart2(abs(Vart2)>pi)=NaN;

bw = pi/6;
bins = [min(Vart1(:)):bw:max(Vart1(:)) 2*pi];
meanfish = NaN(length(bins),length(unique(FishID)));
for b = 1 : length(bins-1)-1
    for fi = unique(FishID)'
        seqs = find(FishID==fi);
        v2 = Vart2(seqs,:);
        v1 = Vart1(seqs,:);
        meanfish(b,fi) = nanmean(v2(v1>bins(b) & v1 <= bins(b+1)))-nanmean(v2(:));
    end
end
%***
figure
p = plot(repmat(bins+bw/2, [length(unique(FishID)), 1])', meanfish);
for i = 1 : length(p)
    p(i).Color = [0.7 0.7 0.7];
    p(i).Marker = 'o';
    p(i).MarkerSize = 8;
    p(i).LineStyle = 'none';
end

fish_in_bin = size(meanfish,2) - sum(isnan(meanfish),2);
hold on
errorbar(bins+bw/2,  nanmean(meanfish, 2), nanstd(meanfish,1,2)./sqrt(fish_in_bin), 'k', 'Linewidth', 1.5)

ylabel('<\delta\theta_n>')
xlabel('\theta_n');
xlim([-0.01 2*pi+0.01])
xticks([0 : pi/2 : 2*pi])
xticklabels({'-\pi','-\pi/2', '0','\pi/2', '\pi'})

ylim([-pi/2 pi/2])
yticks([-pi/2 : pi/2 : pi/2])
yticklabels({'-\pi/2', '0','\pi/2'})


ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';

%% <dX> VS I per fish
Vart1 = Lum(:, 1:end-1);

xl = 'I_n';
Vart2 = dX_corr(:, 1:end); % dX corrected from individual bias
Vart2(abs(Vart2)>pi)=NaN;
b = 9;
yl = '<\delta\theta_n>';

[binvals, elts_per_bin, dXdII] = BinsWithEqualNbofElements(Vart1, Vart2, b, b+3);

%***
f = figure;
xlimits = [0 0.23];
area(xlimits, [sigma_error/2, sigma_error/2], -sigma_error/2,...
     'FaceColor', colour(4,:), 'FaceAlpha', 0.3, 'LineStyle', 'none', 'ShowBaseLine', 'off',...
     'DisplayName', '\sigma_{\theta measurement error}')
hold on
errorbar(binvals,  mean(dXdII, 2), std(dXdII,1,2)/sqrt(elts_per_bin-1),...
    'Color', colour(2,:), 'Linewidth', 1.5, 'DisplayName', 'mean on all temporal experiments')
xlabel(xl)
ylabel(yl)
xlim(xlimits)
ylim([-0.1 0.1])

ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
legend