% Enuc analysis

%%
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Figures/';
%%
Enuc.pooled.Load

%%
subplot(3,2,1)
histogram(dX)

subplot(3,2,2)
histogram(IBI)
text(20,500, {['mean : ' num2str(nanmean(IBI(:)))  ],...
    ['median : ' num2str(nanmedian(IBI(:)))  ]});
title('IBI')

subplot(3,2,3)
histogram(Dist)
text(150,300, {['mean : ' num2str(nanmean(Dist(:)))  ],...
    ['median : ' num2str(nanmedian(Dist(:)))  ]});
title('distance/bout (px)')

%% <dX^2>n &vs dI/I(n-1)
Vart1_b = dLum(:, 1:end-1)./( (Luminosity(:, 1:end-2)+ Luminosity(:, 2:end-1))./2 ) ;
xl = '\deltaI/I_{n-1}';
Vart2 = dX(:, 2:end).^2;
b = 12;
yl = '<\delta\theta^2>_n';

[binvals, elts_per_bin, v2bin, ~, binedges] = BinsWithEqualNbofElements(Vart1_b, Vart2, b, b+10);

iszero = find(binvals==0);
v2bin_onezerobin = v2bin;
v2bin_onezerobin(iszero(1),:)= mean(v2bin(iszero,:));
v2bin_onezerobin(iszero(2:end),:) = [];
binvals(iszero(2:end)) = [];
binedges(iszero(2:end)) = [];

%***
f = figure;
errorbar(binvals,  mean(v2bin_onezerobin, 2), std(v2bin_onezerobin,1,2)/sqrt(elts_per_bin), std(v2bin_onezerobin,1,2)/sqrt(elts_per_bin),...
    diff(binedges)/4, diff(binedges)/4, 'k', 'Linewidth', 1.5)
hold on
plot([binvals(1)*1.8, binvals(end)*1.8], [0.18 0.18], 'k--')
xlim([min(binvals)*1.8 max(binvals)*1.8])
xticks([-1:0.5:1.5])
xlabel(xl)
ylabel(yl)
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%--- stat test ---
h = NaN(size(v2bin_onezerobin,1));
p = NaN(size(v2bin_onezerobin,1));
for i = 1 : size(v2bin_onezerobin,1)
    for j = 1 : size(v2bin_onezerobin,1)
        [h(i,j), p(i,j)] = ttest2(v2bin_onezerobin(i,:),v2bin_onezerobin(j,:));
    end
end

%***
figure
heatmap(binvals, binvals, round(p,2))

colormap('gray')

%% BUILT-IN double gaussian fit on dX
boutsperseq = size(X,2)-sum(isnan(X),2);
mednbouts = median(boutsperseq);

% variables
Vart1 = dLum(:, 1:end-1)./( (Luminosity(:, 1:end-2)+ Luminosity(:, 2:end-1))./2 ) ;
Vart2 = dX(:, 2:end);
dX(dX==0)=NaN;
[binvals, elts_per_bin, v2bin_pos] = BinsWithEqualNbofElements(Vart1, Vart2, 3, 7);

f1 = NaN(size(v2bin_pos,1), 3);
f2  = NaN(size(v2bin_pos,1), 3);
nf = NaN(size(v2bin_pos,1), 1);
nt = NaN(size(v2bin_pos,1), 1);
ci1low = NaN(size(v2bin_pos,1), 3);
ci1high = NaN(size(v2bin_pos,1), 3);
ci2low = NaN(size(v2bin_pos,1), 3);
ci2high = NaN(size(v2bin_pos,1), 3);
for i = 1 : size(v2bin_pos,1)
    c = range(v2bin_pos(i,:));
    if c == 0
        c =1;
    end
    [dx_histogram, x_histogram] = hist(v2bin_pos(i,:), 2*round(sqrt(elts_per_bin)));
    f = fit(x_histogram.',dx_histogram.','gauss2');
    forward = @(x) f.a1*exp(-((x-f.b1)/f.c1).^2);
    side = @(x) f.a2*exp(-((x-f.b2)/f.c2).^2);
        figure;
        plot(x_histogram, dx_histogram)
        hold on
        plot(x_histogram, f.a1*exp(-((x_histogram-f.b1)/f.c1).^2) + f.a2*exp(-((x_histogram-f.b2)/f.c2).^2))
        xlim([-pi pi])
    f1(i, 1) = f.a1;
    f1(i, 2) = f.b1;
    f1(i, 3) = f.c1;
    f2(i, 1) = f.a2;
    f2(i, 2) = f.b2;
    f2(i, 3) = f.c2;
    nf(i) = sum(forward(x_histogram));
    nt(i) = sum(side(x_histogram));
    ci = confint(f);
    ci1low(i,1:3) = ci(1,1:3);
    ci1high(i,1:3) = ci(2,1:3);
    
    f2(i, 1) = f.a2;
    f2(i, 2) = f.b2;
    f2(i, 3) = f.c2;
    ci2low(i,1:3) = ci(1,4:6);
    ci2high(i,1:3) = ci(2,4:6);
end

%% ***
%  (2) <dX2>
%***
fig1 = figure;
yyaxis left
errorbar(binvals/max(abs(binvals)), f2(:,3).^2, f2(:,3).^2-ci2low(:,3).^2, -f2(:,3).^2+ci2high(:,3).^2,...
    'k-', 'Linewidth', 1.5)
hold on
plot([binvals(1)*1.05 binvals(end)*1.05], [0.5476 0.5476], 'k--')
ylabel('<\delta\theta^2>_t_u_r_n')
ylim([0 5])
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'normal';
ax.YColor = [0 0 0];
yyaxis right
errorbar(binvals/max(abs(binvals)), f1(:,3).^2, f1(:,3).^2-ci1low(:,3).^2, -f1(:,3).^2+ci1high(:,3).^2,...
    'Color', [0.3 0.6 0.7], 'Linewidth', 1.5)
hold on
plot([binvals(1)*1.1 binvals(end)*1.1], [0.0121 0.0121], 'Color', [0.3 0.6 0.7])
ylabel('<\delta\theta^2>_s_c_o_o_t')
ylim([0 0.5])
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'bold';
ax.YColor = [0.3 0.6 0.7];

xlim([binvals(1)*0.8 1])
xlabel('relative contrast (I_L-I_R)/I_m_a_x')
title({'Mean bout variance',...
    'from double gaussian fit'}')

%  (3) P(turn)
%***
fig2 = figure;
plot(binvals/max(abs(binvals)), f2(:,1).*f2(:,3)./(f1(:,1).*f1(:,3)+f2(:,1).*f2(:,3)),...
    'k-sq', 'Markersize', 3, 'MarkerFaceColor', 'k')
hold on
plot([binvals(1) 1], [0.49 0.49], 'k--')
xlim([-0.8 1.2])
ylim([0 1])
ylabel('P(turn)')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'bold';
title('')
