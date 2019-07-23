

colour = colour_palette(0,4);

%% BUILT-IN double gaussian fit on dX

% variables
Var1 = dLum(:, 1:mednbouts)./( (Lum(:, 1:mednbouts)+ Lum(:, 2:mednbouts+1))./2 ) ;
Var2 = dX(:, 2:mednbouts+1);
dX(dX==0)=NaN;
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Var1, Var2, 21, 24);

f1 = NaN(size(v2bin,1), 3);
f2  = NaN(size(v2bin,1), 3);
nf = NaN(size(v2bin,1), 1);
nt = NaN(size(v2bin,1), 1);
ci1low = NaN(size(v2bin,1), 3);
ci1high = NaN(size(v2bin,1), 3);
ci2low = NaN(size(v2bin,1), 3);
ci2high = NaN(size(v2bin,1), 3);
%***
figure
hold on
for i = 1 : size(v2bin,1)
    c = range(v2bin(i,:));
    if c == 0
        c =1;
    end
    [dx_histogram, x_histogram] = hist(v2bin(i,:), c/2*round(sqrt(elts_per_bin)));
    f = fit(x_histogram.',dx_histogram.','gauss2');
    forward = @(x) f.a1*exp(-((x-f.b1)/f.c1).^2);
    side = @(x) f.a2*exp(-((x-f.b2)/f.c2).^2);
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
    ci = confint(f, 0.99);
    ci1low(i,1:3) = ci(1,1:3);
    ci1high(i,1:3) = ci(2,1:3);
    
    f2(i, 1) = f.a2;
    f2(i, 2) = f.b2;
    f2(i, 3) = f.c2;
    ci2low(i,1:3) = ci(1,4:6);
    ci2high(i,1:3) = ci(2,4:6);
end

%  (2) <dX2>
%***
fig1 = figure;
errorbar(binvals/max(abs(binvals)), f2(:,3), f2(:,3)-ci2low(:,3), -f2(:,3)+ci2high(:,3),...
    'k-', 'Linewidth', 1.5, 'DisplayName', '\sigma_{turn} \pm ci99%')
hold on
plot([binvals(1)*1.05 binvals(end)*1.05], [0.69 0.69],...
    'k--', 'DisplayName', '\sigma_{turn} neutral')
errorbar(binvals/max(abs(binvals)), f1(:,3), f1(:,3)-ci1low(:,3), -f1(:,3)+ci1high(:,3),...
    'Color', colour(4,:), 'Linewidth', 1.5, 'DisplayName', '\sigma_{scoot} \pm ci99%')
plot([binvals(1)*1.1 binvals(end)*1.1], [0.12 0.12],...
    '--', 'Color', colour(2,:), 'DisplayName', '\sigma_{scoot} neutral')

ylabel('<\sigma>')
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

xlim([binvals(1)*1.1 binvals(end)*1.1])
xlabel('relative contrast (I_L-I_R)/I_m_a_x')
title({'Mean bout variance',...
    'from double gaussian fit'}')

legend

%%  (3) P(turn)
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

%% custom gaussian fit
% variables
Var1 = dLum(:, 1:mednbouts)./( (Lum(:, 1:mednbouts)+ Lum(:, 2:mednbouts+1))./2 ) ;
Var2 = dX(:, 2:mednbouts+1);

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Var1, Var2, 12, 18);

wturn = nan(1, size(v2bin,1));
wfor = nan(1, size(v2bin,1));
pturn = nan(1, size(v2bin,1));
pfor = nan(1, size(v2bin,1));
sep = nan(size(v2bin,1),1);
sewt = nan(size(v2bin,1),1);
sewf = nan(size(v2bin,1),1);

for i = 1 : size(v2bin,1)
    selection = v2bin(i,:);
    selection(abs(selection)>pi) = NaN;
    [wturn(i), wfor(i), pturn(i), pfor(i), sep(i), sewt(i), sewf(i)] = Gauss2custom(selection );
end
    
%  (2) <dX2>
%***
fig = figure;
errorbar(binvals/max(abs(binvals)), wturn, sewt, ...
    'sq-', 'Color', colour(2,:), 'Linewidth', 1.5)
hold on
errorbar(binvals/max(abs(binvals)), wfor, sewf,...
    'sq-', 'Color', colour(4,:), 'Linewidth', 1.5)
ylim([0 1.5])
ylabel('<\delta\theta^2>')
xlim([binvals(1)*1.1 binvals(end)*1.1])
xlabel('relative contrast dI/I')
title({'Mean bout variance',...
    'from custom double gaussian fit'}')

legend('<\delta\theta^2>_{turn}', '<\delta\theta^2>_{scoot}')

ax=gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';

%  (3) P(turn)
%***
figure
plot([-0.85 1.02], [0.4089 0.4089], '--', 'LineWidth', 1.5, 'Color', colour(5,:))
hold on
errorbar(binvals/max(abs(binvals)), pturn, sep,...
    'k-sq', 'Markersize', 6, 'MarkerFaceColor', colour(2,:), 'LineWidth', 1.5)
legend('Neutral', 'Temporal stimulus')
ylim([0 1])
ylabel('P(turn)')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'bold';
title('')
ax=gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';

%% --- fit the distribution of dX with a CUSTOM double-gaussian with constrained mean and sigma ---
colour = colour_palette(0,4);
colour_spont = colour_palette(0,1);

Var1 = dLum(:, 1:mednbouts)./( (Lum(:, 1:mednbouts)+ Lum(:, 2:mednbouts+1))./2 ) ;
%Var1 = dLum(:, 1:end-1)./( (Lum(:, 1:end-2)+ Lum(:, 2:end-1))/2 ) ;
Var1(Var1 == 0)=NaN;

Var2 = dX(:, 2:mednbouts+1);

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Var1, Var2, 12, 18);
mV2 = nanmean(v2bin,2);
a = nan(size(v2bin,1),1);
mut = nan(size(v2bin,1),1);
sigmat = nan(size(v2bin,1),1);
aconfint = nan(size(v2bin,1),2);
muconfint = nan(size(v2bin,1),2);
sigmaconfint = nan(size(v2bin,1),2);

figure
muf = 0;
sigmaf = 0.12;
for i = 1 : size(v2bin, 1)
    v = (v2bin(i,:)); 
    v = v(:);
    v(isnan(v)) = [];
    N = length(v);
    binwidth = 2*iqr(v)/(N^(1/3));
    bins = min(v) : binwidth : max(v);
    [dx_histogram, x_histogram]  = histcounts(v,bins);
    x_histogram = x_histogram(1:end-1) + binwidth/2;
    dx_histogram = dx_histogram/(sum(dx_histogram)*binwidth);
    myfittype = fittype('sum_gauss_constr(x, a, mut, sigmat, muf, sigmaf)',...
        'coefficients', {'a','sigmat','mut'}, 'problem', {'muf', 'sigmaf'});
    myfit = fit(x_histogram',dx_histogram',myfittype, 'startpoint', [0.5 0.7 0], 'problem', {muf, sigmaf});
    
    plot(myfit, x_histogram, dx_histogram)
    %waitforbuttonpress
    a(i) = myfit.a;
    mut(i) = myfit.mut;
    sigmat(i) = myfit.sigmat;
    sep = confint(myfit,0.95);
    aconfint(i,:) = sep(:,1)';
    sigmaconfint(i,:) = sep(:,2)';
    muconfint(i,:) = sep(:,3)';
end

figure
subplot(3,1,1)
errorbar(binvals, mut, mut-muconfint(:,1), muconfint(:,2)-mut,...
    'LineWidth', 1.5, 'Color', colour(3,:), 'DisplayName', '\mu_{turn} estimate')
hold on
plot(binvals, muf*ones(length(binvals),1),...
    '--', 'Linewidth', 1.5, 'Color', colour(4,:), 'DisplayName', '\mu_{forward} fixed')
xlabel('dI/I_{n-1}')
title('\mu')
ylim auto
legend
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

subplot(3,1,2)
plot(binvals, 0.6*ones(length(binvals),1),...
    '-.', 'Linewidth', 1, 'Color', colour_spont(2,:), 'DisplayName', 'spontaneous value')
hold on
errorbar(binvals, sigmat, sigmat-sigmaconfint(:,1), sigmaconfint(:,2)-sigmat,...
    'LineWidth', 1.5, 'Color', colour(3,:), 'DisplayName', '\sigma_{turn} estimate')
hold on
plot(binvals, sigmaf*ones(length(binvals),1),...
    '--', 'Linewidth', 1.5, 'Color', colour(4,:), 'DisplayName', '\sigma_{forward} fixed')
xlabel('dI/I_{n-1}')
title('\sigma')
legend
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

subplot(3,1,3)
plot(binvals, 0.42*ones(length(binvals),1),...
    '-.', 'Linewidth', 1, 'Color', colour_spont(2,:), 'DisplayName', 'spontaneous value')
hold on
errorbar(binvals, a, a-aconfint(:,1),aconfint(:,2)-a,...
    'LineWidth', 1.5, 'Color', colour(4,:), 'HandleVisibility', 'off')
xlabel('dI/I_{n-1}')
title('p_{turn} estimate')
ylim([0 1])
legend
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
