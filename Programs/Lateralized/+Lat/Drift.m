
colour = colour_palette(0,3);
%% Drift of all data
medboutsperseq = median(size(XLat,2) - sum(isnan(XLat),2));

relC = DIlr(:,1:medboutsperseq);
Var1 = relC;
Var2 = (dXLssbiais(:,1:medboutsperseq));
b=7;

% --- bins with equal number of elements
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Var1, Var2, b, b+3);
mV2 = nanmean(v2bin,2);
mV2sq = nanmean(v2bin.^2,2);

% --- linear fit of bias ---
binvals = binvals/max(abs(DIlr(:)));
interval = 1:size(binvals);
linfittype = fittype('a*x','coefficients', {'a'});
myfit = fit(binvals(interval),mV2(interval),linfittype);
Acoeff = myfit.a;

%***
figure;

errorbar(binvals, mV2, std(v2bin,1,2)/sqrt(elts_per_bin),...
    'Linewidth', 2, 'DisplayName', 'Data', 'Color', [0.2 0.2 0.2])
hold on
plot(binvals(interval), Acoeff*binvals(interval),...
    'DisplayName', ['Linear fit : a = ' num2str(Acoeff)], 'Color', colour(1,:), 'Linewidth', 1.5);
ylim([-0.25 0.25])
title(['bias until bout ' num2str(medboutsperseq)])
ylabel('<\delta\theta>')
xlabel('relative contrast')
legend
ax = gca;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';

%***
figure
subplot(1,2,1)
errorbar(binvals, mV2sq, std(v2bin.^2,1,2)/sqrt(elts_per_bin-1),...
    'Linewidth', 2, 'Color', [0.2 0.2 0.2])
ylabel('<\delta\theta^2>')
title('variance')
xlabel('relative contrast')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

subplot(1,2,2)
errorbar(binvals, nanmean(abs(v2bin),2), std(abs(v2bin).^2,1,2)/sqrt(elts_per_bin-1),...
    'Linewidth', 2, 'Color', [0.2 0.2 0.2])
ylabel('<|\delta\theta|>')
title('abs')
xlabel('relative contrast')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

% --- significance test ---
signif_val = 1e-12;
[pm, pv] = mean_var_testing(Var1, Var2, b, signif_val)

%% Drift of all data : ABSOLUTE contrast
relC = DIlr(:,1:17);
absC = abs(relC);
Var1 = absC;
Var2 = (dXLssbiais(:,1:medboutsperseq));
Var2(relC<0) = - Var2(relC<0);
b=6;

% --- bins with equal number of elements
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Var1, Var2, b, b+3);
mV2 = nanmean(v2bin,2);
mV2sq = nanmean(v2bin.^2,2);

% --- linear fit of bias ---
binvals = binvals/max(abs(DIlr(:)));
interval = 1:size(binvals);
linfittype = fittype('a*x','coefficients', {'a'});
myfit = fit(binvals(interval),mV2(interval),linfittype);
Acoeff = myfit.a;

%***
figure;

errorbar(binvals, mV2, std(v2bin,1,2)/sqrt(elts_per_bin),...
    'Linewidth', 2, 'DisplayName', 'Data', 'Color', [0.2 0.2 0.2])
hold on
plot(binvals(interval), Acoeff*binvals(interval),...
    'DisplayName', ['Linear fit : a = ' num2str(Acoeff)], 'Color', colour(1,:), 'Linewidth', 1.5);
ylim([-0.1 0.3])
title('bias')
ylabel('<\delta\theta>')
xlabel('absolute contrast')
legend
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%***
figure
subplot(1,2,1)
errorbar(binvals, mV2sq, std(v2bin.^2,1,2)/sqrt(elts_per_bin-1),...
    'Linewidth', 2, 'Color', [0.2 0.2 0.2])
ylabel('<\delta\theta^2>')
title('variance')
xlabel('absolute contrast')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

subplot(1,2,2)
errorbar(binvals, nanmean(abs(v2bin),2), std(abs(v2bin).^2,1,2)/sqrt(elts_per_bin-1),...
    'Linewidth', 2, 'Color', [0.2 0.2 0.2])
ylabel('<|\delta\theta|>')
title('abs')
xlabel('absolute contrast')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%% Drift on data per fish
relC = DIlr(:, 1:end-1);
Var1 = relC;
Var2 = (dXLssbiais(:,1:medboutsperseq));
Var2 = (dXLssbiais);

b=12;

fish = unique(FishID);
alpha = NaN(length(fish),1);
alphaci = NaN(length(fish),2);
medboutsperfish = NaN(length(fish),1);
minboutsperfish = NaN(length(fish),1);
maxboutsperfish = NaN(length(fish),1);

figure;
hold on
for i = 1:length(fish)
    seqafish = (FishID == i);
    var1 = Var1(seqafish,:);
    var2 = Var2(seqafish,:);
    if length(var1(:)) - sum(isnan(var1(:))) < b+2
        continue
    end
    [binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, b, b+2);
    mV2 = nanmean(v2bin,2);
    
    plot(binvals, mV2)
    linfittype = fittype('a*x','coefficients', {'a'});
    myfit = fit(binvals, mV2, linfittype);
    alpha(i) = myfit.a;
    alphaci(i,:) = confint(myfit)';
    medboutsperfish(i) = median(size(var1,2) - sum(isnan(var1),2));
    minboutsperfish(i) = min(size(var1,2) - sum(isnan(var1),2));
    maxboutsperfish(i) = max(size(var1,2) - sum(isnan(var1),2));
end

%***
figure;
plot(fish, alpha, 'o')
legend
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

figure;
errorbar(medboutsperfish, alpha, alpha-alphaci(:,1), alphaci(:,2)-alpha, minboutsperfish, maxboutsperfish,...
    'o', 'LineWidth', 1.5, 'Color', 'k', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
ax = gca;
ax.XScale = 'log';
ax.FontSize = 16;
ax.FontName = 'Times New Roman';
xlabel('median number of bouts per fish')

figure
plot(alpha, Rprojperfish, 'sq')
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Times New Roman';


%% Drift of all data : THETA
ThetaSel = XLat(:,1:medboutsperseq);
ThetaSel = wrapTo2Pi(ThetaSel)-pi;
Var1 = ThetaSel;
Var2 = (dXLssbiais(:,1:medboutsperseq));
b=12;

% --- bins with equal number of elements
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Var1, Var2, b, b+3);
mV2 = nanmean(v2bin,2);
mV2sq = nanmean(v2bin.^2,2);

% --- linear fit of bias ---
find_peak = find(binvals>=0, 1);
interval1 = 1 : find_peak;
interval2 = find_peak - 1 : size(binvals);
linfittype = fittype('a*x','coefficients', {'a'});
myfit1 = polyfit(binvals(interval1), mV2(interval1), 1);
myfit2 = polyfit(binvals(interval2), mV2(interval2), 1);
acoeff1 = myfit1(1);
acoeff2 = myfit2(1);

%***
figure;
errorbar(binvals, mV2, std(v2bin,1,2)/sqrt(elts_per_bin),...
    'Linewidth', 2, 'DisplayName', 'Data', 'Color', [0.2 0.2 0.2])
hold on
plot(binvals(interval1), myfit1(2)+acoeff1*binvals(interval1),...
    'DisplayName', ['Linear fit : a = ' num2str(acoeff1)], 'Color', colour(1,:), 'Linewidth', 1.5);
plot(binvals(interval2), myfit2(2)+acoeff2*binvals(interval2),...
    'DisplayName', ['Linear fit : a = ' num2str(acoeff2)], 'Color', colour(1,:), 'Linewidth', 1.5);
title('bias')
ylabel('<\delta\theta>')
xlabel('\theta')
legend
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%***
figure
subplot(1,2,1)
errorbar(binvals, mV2sq, std(v2bin.^2,1,2)/sqrt(elts_per_bin-1),...
    'Linewidth', 2, 'Color', [0.2 0.2 0.2])
ylabel('<\delta\theta^2>')
title('variance')
xlabel('absolute contrast')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

subplot(1,2,2)
errorbar(binvals, nanmean(abs(v2bin),2), std(abs(v2bin).^2,1,2)/sqrt(elts_per_bin-1),...
    'Linewidth', 2, 'Color', [0.2 0.2 0.2])
ylabel('<|\delta\theta|>')
title('abs')
xlabel('absolute contrast')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%% --- drift in time (A)
tbinning = 1:50;
win = 1;

Var1 = DIlr(:,1:end-1)/(max(abs(DIlr(:))));
Var2 = dXl;

acoeff = NaN(size(tbinning));
for i = 1 : length(tbinning)-1
    vart1 = Var1(:,tbinning(i):tbinning(i)+win);
    vart2 = Var2(:,tbinning(i):tbinning(i)+win);
    b=7;
    
    % --- bins with equal number of elements
    [binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(vart1, vart2, b, b+10);
    mV2 = nanmean(v2bin,2);
    plot(binvals, mV2)
    hold on
    % --- linear fit of bias ---
    interval = 1:size(binvals);
    linfittype = fittype('a*x','coefficients', {'a'});
    myfit = fit(binvals(interval),mV2(interval),linfittype);
    acoeff(i) = myfit.a;
end
figure;
plot(acoeff, 'LineWidth', 1.5)
hold on
plot(smooth(acoeff,5), 'LineWidth', 1.5)
xlabel('bout #')
ylabel('\alpha')
title(['bias coefficient on ' num2str(win) ' bouts window (smoothed over 5 bouts)'])
ax = gca;
ax.XScale = 'lin';
ax.FontSize = 16;
ax.FontName = 'Times New Roman';

%% --- fit the distribution of dX with a double-gaussian ---
Var1 = DIlr(:, 1:17)/(max(abs(DIlr(:))));
Var2 = dXLssbiais(:,1:17);

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Var1, Var2, 7, 9);
mV2 = nanmean(v2bin,2);
mV2sq = nanvar(v2bin,1,2);

f1 = NaN(size(v2bin,1), 3);
f2  = NaN(size(v2bin,1), 3);
ci1 = NaN(size(v2bin,1), 3, 2);
ci2 = NaN(size(v2bin,1), 3, 2);
std1 = NaN(size(v2bin,1), 3);
std2 = NaN(size(v2bin,1), 3);
%***
figure
for i = 1 : size(v2bin, 1)
    [dx_histogram, x_histogram] = hist(v2bin(i,:), 2*round(sqrt(elts_per_bin)),...
        'Normalization', 'probability');
    [f, lol, lil] = fit(x_histogram.',dx_histogram.','gauss2');
    forward = @(x) f.a1*exp(-((x-f.b1)/f.c1).^2);
    side = @(x) f.a2*exp(-((x-f.b2)/f.c2).^2);
    
    plot(x_histogram, dx_histogram)
    hold on
    plot(x_histogram, f.a1*exp(-((x_histogram-f.b1)/f.c1).^2) + f.a2*exp(-((x_histogram-f.b2)/f.c2).^2))
    %waitforbuttonpress
    
    CI = confint(f, 0.99);
    % forward
    f1(i, 1) = f.a1;
    f1(i, 2) = f.b1;
    f1(i, 3) = f.c1;
    ci1(i,1:3,:) = CI(:,1:3)';
    
    % turn
    f2(i, 1) = f.a2;
    f2(i, 2) = f.b2;
    f2(i, 3) = f.c2;
    ci2(i,1:3,:) = CI(:,4:6)';
    
    CI = confint(f, 0.68);
    std1(i,:) = (CI(2,1:3)-CI(1,1:3))/2;
    std2(i,:) = (CI(2,4:6)-CI(1,4:6))/2;
end

%***
fig = figure;
vc = [255/255,127/255,80/255];
yyaxis left
errorbar(binvals, mV2, std(v2bin,1,2)/sqrt(elts_per_bin), 'Linewidth', 2)
ylabel('dX')

yyaxis right
errorbar(binvals, mV2sq, std(v2bin.^2,1,2)/sqrt(elts_per_bin), 'b', 'Linewidth', 1.5, 'Color', vc)
ylabel('dX^2')
ylim([0 1])
grid on
ax=gca;
ax.FontSize = 14;

%***
fig = figure;
subplot(2,1,1)
errorbar(binvals, f1(:,2), ci1(:,2,2)-f1(:,2), f1(:,2)-ci1(:,2,1),...
    '--', 'LineWidth', 1.5, 'Color', colour(4,:), 'DisplayName', 'forward distribution (99%ci)')
hold on
errorbar(binvals, f2(:,2), ci2(:,2,2)-f2(:,2), f2(:,2)-ci2(:,2,1),...
    'LineWidth', 1.5, 'Color', colour(2,:), 'DisplayName', 'turn distribution(99%ci)')
ylabel('\mu_{d\theta}')
xlabel('C_{rel}')
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 14;
legend

subplot(2,1,2)
errorbar(binvals, f1(:,3), ci1(:,3,2)-f1(:,3), f1(:,3)-ci1(:,3,1),...
    '--', 'LineWidth', 1.5, 'Color', colour(4,:), 'DisplayName', 'forward distribution (99%ci)')
hold on
errorbar(binvals, f2(:,3), ci2(:,3,2)-f2(:,3), f2(:,3)-ci2(:,3,1),...
    'LineWidth', 1.5, 'Color', colour(2,:), 'DisplayName', 'turn distribution (99%ci)')
xlabel('C_{rel}')
ylabel('\sigma_{d\theta}')
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 14;

pturn = f2(:,1).*f2(:,3) ./ ( f1(:,1).*f1(:,3)+f2(:,1).*f2(:,3));
stdpturn = pturn.*sqrt( (std2(:,1)./f2(:,1)).^2 + (std2(:,3)./f2(:,3)).^2 + (std1(:,1)./f1(:,1)).^2 + (std2(:,3)./f2(:,3)).^2 );

%***
fig2 = figure;
errorbar(binvals,  pturn, stdpturn,...
    'k-sq', 'Markersize', 6, 'MarkerFaceColor', 'k')
ylim([0 1])
xlabel('IL-IR')
ylabel('P(turn)')
title('turns/forward')


%% --- fit the distribution of dX with a CUSTOM double-gaussian ---
Var1 = DIlr(:,1:17)/(max(abs(DIlr(:))));
Var2 = dXLssbiais(:,1:17);

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Var1, Var2, 7, 9);

mV2 = nanmean(v2bin,2);
sV2 = nanstd(v2bin,1,2);
maV2 = nanmean(abs(v2bin),2);
saV2 = nanstd(abs(v2bin),1,2);

a = nan(size(v2bin,1),1);
mut = nan(size(v2bin,1),1);
sigmat = nan(size(v2bin,1),1);
muf = nan(size(v2bin,1),1);
sigmaf = nan(size(v2bin,1),1);
aconfint = nan(size(v2bin,1),2);
aconfint68 = nan(size(v2bin,1),2);
muconfint = nan(size(v2bin,1),2);
sigmaconfint = nan(size(v2bin,1),2);
figure
for i = 1 : size(v2bin, 1)
    v = (v2bin(i,:)); 
    v = v(:);
    v(isnan(v)) = [];
    N = length(v);
    binwidth = 3*iqr(v)/(N^(1/3));
    bins = min(v) : binwidth : max(v);
    [dx_histogram, x_histogram]  = histcounts(v,bins);
    x_histogram = x_histogram(1:end-1) + binwidth/2;
    dx_histogram = dx_histogram/(sum(dx_histogram)*binwidth);
    muf = 0;
    sigmaf = 0.12;
    myfittype = fittype('sum_gauss_constr(x, a, mut, sigmat, muf, sigmaf)',...
        'coefficients', {'a','mut','sigmat'}, 'problem', {'muf', 'sigmaf'}); %
    myfit = fit(x_histogram',dx_histogram',myfittype, 'startpoint', [0.5 0 0.7], 'Upper', [1 Inf 0.9], 'problem', {muf, sigmaf});%
    
    plot(myfit, x_histogram, dx_histogram)
    a(i) = myfit.a;
    mut(i) = myfit.mut;
    sigmat(i) = myfit.sigmat;
%     muf(i) = myfit.muf;
%     sigmaf(i) = myfit.sigmaf;
    ci = confint(myfit,0.99);
    aconfint(i,:) = ci(:,1)';
    muconfint(i,:) = ci(:,2)';
    sigmaconfint(i,:) = ci(:,3)';
    
    ci = confint(myfit,0.68);
    aconfint68(i,:) = ci(:,1)';
end
%plot(binvals, a.*sigmat./(a.*sigma1+a2.*sigma2))
figure
subplot(3,1,1)
errorbar(binvals, mut, mut-muconfint(:,1), muconfint(:,2)-mut,...
    'LineWidth', 1.5, 'Color', colour(3,:), 'DisplayName', '\mu_{turn} estimate')
hold on
plot(binvals, muf*ones(length(binvals),1),...
    '--', 'Linewidth', 1.5, 'Color', colour(4,:), 'DisplayName', '\mu_{forward} fixed')
xlabel('C_{rel}')
title('\mu from custom fit')
ylim auto
legend

subplot(3,1,2)
errorbar(binvals, sigmat, sigmat-sigmaconfint(:,1), sigmaconfint(:,2)-sigmat,...
    'LineWidth', 1.5, 'Color', colour(3,:), 'DisplayName', '\sigma_{turn} estimate')
hold on
plot(binvals, sigmaf*ones(length(binvals),1),...
    '--', 'Linewidth', 1.5, 'Color', colour(4,:), 'DisplayName', '\sigma_{forward} fixed')
xlabel('C_{rel}')
title('\sigma from custom fit')
legend

subplot(3,1,3)
errorbar(binvals, a, a-aconfint(:,1),aconfint(:,2)-a,...
    'LineWidth', 1.5, 'Color', colour(4,:))
xlabel('C_{rel}')
title('p_{turn} estimate from custom fit')
ylim([0 1])

ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 20;

% percentage of turns in one direction
w = mV2./a./maV2;
dw = w.* sqrt ( (sV2./mV2).^2 +  ( (a-aconfint68(:,1))/2./a ).^2 + (saV2./maV2).^2 );

%***
figure
errorbar(binvals, w, dw/sqrt(elts_per_bin),...
    'LineWidth', 1.5, 'Color', colour(3,:))
xlabel('C_{rel}')
ylabel('<\delta\theta>/(p_{turn}<|\delta\theta|>)')
title('proportion of biased turns')
ylim([-1.2 1.2])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 20;

%% --- fit the distribution of dX with a CUSTOM double-gaussian (2) : ABSOLUTE CONTRAST---
Var1 = DIlr(:, 1:17)/(max(abs(DIlr(:))));
Var2 = dXLssbiais(:,1:17);
Var2(Var1<0) = -Var2(Var1<0);
Var1 = abs(Var1);

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Var1, Var2, 5, 8);

mV2 = nanmean(v2bin,2);

a = nan(size(v2bin,1),1);
mut = nan(size(v2bin,1),1);
sigmat = nan(size(v2bin,1),1);
muf = nan(size(v2bin,1),1);
sigmaf = nan(size(v2bin,1),1);
aconfint = nan(size(v2bin,1),2);
muconfint = nan(size(v2bin,1),2);
sigmaconfint = nan(size(v2bin,1),2);
%***
figure
for i = 1 : size(v2bin, 1)
    v = (v2bin(i,:));
    v = v(:);
    v(isnan(v)) = [];
    N = length(v);
    binwidth = 3*iqr(v)/(N^(1/3));
    bins = min(v) : binwidth : max(v);
    [dx_histogram, x_histogram]  = histcounts(v,bins);
    x_histogram = x_histogram(1:end-1) + binwidth/2;
    dx_histogram = dx_histogram/(sum(dx_histogram)*binwidth);
    muf = 0;
    sigmaf = 0.12;
    myfittype = fittype('sum_gauss_constr(x, a, mut, sigmat, muf, sigmaf)',...
        'coefficients', {'a','mut', 'sigmat'}, 'problem', {'muf', 'sigmaf'}); %
    myfit = fit(x_histogram',dx_histogram',myfittype, 'startpoint', [0.5 0 0.7], 'problem', {muf, sigmaf});%
    
    plot(myfit, x_histogram, dx_histogram)
    %waitforbuttonpress
    a(i) = myfit.a;
    mut(i) = myfit.mut;
    sigmat(i) = myfit.sigmat;
%     muf(i) = myfit.muf;
%     sigmaf(i) = myfit.sigmaf;
    ci = confint(myfit,0.99);
    aconfint(i,:) = ci(:,1)';
    muconfint(i,:) = ci(:,2)';
    sigmaconfint(i,:) = ci(:,3)';
end
%plot(binvals, a.*sigmat./(a.*sigma1+a2.*sigma2))
figure
subplot(3,1,1)
plot(binvals, muf*ones(length(binvals),1),...
    '--', 'Linewidth', 1.5, 'Color', colour(4,:), 'DisplayName', '\mu_{forward} fixed')
hold on
errorbar(abs(binvals), mut, mut-muconfint(:,1), muconfint(:,2)-mut,...
    'LineWidth', 1.5, 'Color', colour(3,:), 'DisplayName', '\mu_{turn} estimate')
xlabel('C_{abs}')
title('\mu')
ylim auto
legend

subplot(3,1,2)
hold on
plot(binvals, sigmaf*ones(length(binvals),1),...
    '--', 'Linewidth', 1.5, 'Color', colour(4,:), 'DisplayName', '\sigma_{forward} fixed')
errorbar(binvals, sigmat, sigmat-sigmaconfint(:,1), sigmaconfint(:,2)-sigmat,...
    'LineWidth', 1.5, 'Color', colour(3,:), 'DisplayName', '\mu_{forward} estimate')
xlabel('C_{rel}')
title('\sigma')
legend

subplot(3,1,3)
errorbar(binvals, a, a-aconfint(:,1), aconfint(:,2)-a,...
    'LineWidth', 1.5, 'Color', colour(2,:))
xlabel('C_{rel}')
title('p_{turn} estimate ci99%')
ylim([0 1])

%***
figure
errorbar(abs(binvals), mut, mut-muconfint(:,1), muconfint(:,2)-mut,...
    'LineWidth', 1.5, 'Color', colour(3,:), 'DisplayName', '\mu_{turn} estimate')
xlabel('C_{abs}')
title('\mu')
ylim auto
legend


%% create drift matrix
Contrasts = DIlr(:,1:23)/max(abs(DIlr(:)));
A = Contrasts*Acoeff;



