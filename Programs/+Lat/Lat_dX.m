%%


%%
bornes = [1 20];
%x = Xfilt(:, bornes(1) : bornes(2));

x = XLat;%(TimeBout<60);

dx = diff(x, 1, 2 );
%dx(abs(dx)>pi) = NaN;

%%
Vart1 = DIlr(:, 1:end-1);
Vart2 = dx;

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, 6, 10);
mV2 = nanmean(v2bin,2);

%***
fig = figure;
hold on
plot([-400 400], [-0.1 0.1], 'Color', [0.8 0.8 0.8])
plot([-400 400], [0 0], 'Color', [0.8 0.8 0.8])
plot([0 0], [-0.1 0.1], 'Color', [0.8 0.8 0.8])
errorbar(binvals, mV2, std(v2bin,1,2)/sqrt(elts_per_bin), 'k', 'Linewidth', 2)
yticks(-0.1:0.025:0.1)
%xticks(-1:0.5:1)
xticks(-400:100:400)
%xlim([-1 1])
xlim([-400 400])
ylim([-0.1 0.1])
%xlabel('Relative contrast (I_L-I_R)')
xlabel('\DeltaI = I_L-I_R (\muW.cm^-^2)')
ylabel('<dX> (rad)')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%% save figure
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/FiguresSerious/';
name = ['lateralized_dX_abscontrast'];
saveas(fig,[figpath name],'fig')
disp('fig saved')

%% GAUSSIAN FITS
%..........................................................................
%% --- built-in gauss2 : all data ---

Vart2 = dx;

% --- gauss 2 on whole dataset ---
N = dx(:);
N(isnan(N))=[];
N = length(N);
[dx_histogram, x_histogram] = hist(Vart2(:), 2*round((N)^(1/2)),...
    'Normalization', 'probability');
f = fit(x_histogram.',dx_histogram.','gauss2');
forward = @(x) f.a1*exp(-((x-f.b1)/f.c1).^2);
side = @(x) f.a2*exp(-((x-f.b2)/f.c2).^2);

%***
figure;
plot(x_histogram, dx_histogram)
hold on
plot(x_histogram, f.a1*exp(-((x_histogram-f.b1)/f.c1).^2) + f.a2*exp(-((x_histogram-f.b2)/f.c2).^2))
xlim([-pi pi])

%% --- built-in gauss2 : on binned variables ---

Vart1 = DIlr(:, 1:end-1);%wrapToPi(XLat(:, 1:end-1));
Vart2 = dx;

% --- bin variables ---
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, 6, 10);
mV2 = nanmean(v2bin,2); % mean
mV2sq = nanmean(v2bin.^2,2); %variance

%***
fig = figure;
vc = [255/255,127/255,80/255];

yyaxis left
errorbar(binvals/max(abs(binvals)), mV2, std(v2bin,1,2)/sqrt(elts_per_bin), 'Linewidth', 2)
ylabel('dX')

yyaxis right
plot(binvals/max(abs(binvals)), mV2sq, 'b', 'Linewidth', 1.5, 'Color', vc)
ylabel('dX^2')
grid on
ax=gca;
ax.FontSize = 14;


% --- gauss fit on bins ---
f1 = NaN(size(v2bin,1), 3);
f2  = NaN(size(v2bin,1), 3);
ci1low = NaN(size(v2bin,1), 3);
ci1high = NaN(size(v2bin,1), 3);
ci2low = NaN(size(v2bin,1), 3);
ci2high = NaN(size(v2bin,1), 3);

for i = 1 : size(v2bin, 1)
    [dx_histogram, x_histogram] = hist(v2bin(i,:), round(range(v2bin(i,:))/7*(elts_per_bin)^(1/2)),...
        'Normalization', 'probability');
    disp( round(range(v2bin(i,:))/2*(elts_per_bin)^(1/2)))
    f = fit(x_histogram.',dx_histogram.','gauss2');
    forward = @(x) f.a1*exp(-((x-f.b1)/f.c1).^2);
    side = @(x) f.a2*exp(-((x-f.b2)/f.c2).^2);
        figure;
        plot(x_histogram, dx_histogram)
        hold on
        plot(x_histogram, f.a1*exp(-((x_histogram-f.b1)/f.c1).^2) + f.a2*exp(-((x_histogram-f.b2)/f.c2).^2))
    f1(i, 1) = f.a1;
    f1(i, 2) = f.b1;
    f1(i, 3) = f.c1;
    ci = confint(f);
    ci1low(i,1:3) = ci(1,1:3);
    ci1high(i,1:3) = ci(2,1:3);
    
    f2(i, 1) = f.a2;
    f2(i, 2) = f.b2;
    f2(i, 3) = f.c2;
    ci2low(i,1:3) = ci(1,4:6);
    ci2high(i,1:3) = ci(2,4:6);
end

%% *** (1) <dX> from gauss2
fig = figure;

fig = figure;
hold on
errorbar(binvals/max(abs(binvals)), f2(:,2), f2(:,2)-ci2low(:,2), -f2(:,2)+ci2high(:,2),...
    'Color', 'k', 'Linewidth', 2)
errorbar(binvals/max(abs(binvals)), f1(:,2), f1(:,2)-ci1low(:,2), -f1(:,2)+ci1high(:,2),...
    'Color', [0.3 0.6 0.7], 'Linewidth', 1.5)
ylabel('<dX>')
legend('turn', 'scoot')
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'bold';
title({'Mean bout bias',...
    'from double gaussian fit'})

%% save figure
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/FiguresSerious/';
name = ['lateralized_mean_fromGfit'];
saveas(fig,[figpath name],'fig')
disp('fig saved')
%% ***
%  (2) <dX2>
fig = figure;
xlabel('relative contrast (I_L-I_R)/I_m_a_x')

yyaxis left
errorbar(binvals/max(abs(binvals)), f2(:,3).^2, f2(:,3).^2-ci2low(:,3).^2, -f2(:,3).^2+ci2high(:,3).^2,...
    'k-', 'Linewidth', 2)
ylabel('<dX^2>_t_u_r_n')
ylim([0 2])
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'bold';
ax.YColor = [0 0 0];

yyaxis right
errorbar(binvals/max(abs(binvals)), f1(:,3).^2, f1(:,3).^2-ci1low(:,3).^2, -f1(:,3).^2+ci1high(:,3).^2,...
    'Color', [0.3 0.6 0.7], 'Linewidth', 1.5)
ylabel('<dX^2>_s_c_o_o_t')
ylim([0 0.1])
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'bold';
ax.YColor = [0.3 0.6 0.7];

title({'Mean bout variance',...
    'from double gaussian fit'}')
%% save figure
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/FiguresSerious/';
name = ['lateralized_var_fromGfit'];
saveas(fig,[figpath name],'fig')
disp('fig saved')
%% ***
%  (3) P(turn)
fig = figure;
plot(binvals/max(abs(binvals)), f2(:,1).*f2(:,3)./(f1(:,1).*f1(:,3)+f2(:,1).*f2(:,3)),...
    'k-sq', 'Markersize', 6, 'MarkerFaceColor', 'k')
ylim([0 1])
ylabel('P(turn)')
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'bold';
title('')

%% save figure
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/FiguresSerious/';
name = ['lateralized_pturn_fromGfit'];
saveas(fig,[figpath name],'fig')
disp('fig saved')

%% --- custom double gaussian fit ---

% ::::::::: INCOMPLETE ::::::::::

f1 = NaN(size(v2bin,1), 3);
f2  = NaN(size(v2bin,1), 3);
for i = 1 : size(v2bin, 1)
    [dx_histogram, x_histogram] = hist(v2bin(i,:), 3*round((elts_per_bin)^(1/2)),...
        'Normalization', 'probability');
[alpha, beta, wturn, wfor, pturn, pfor] = Gauss2custom_simple(data);
end

