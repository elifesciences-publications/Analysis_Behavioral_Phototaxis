
%% Variance analysis of temporal phototaxis data
% dXsq versus I ; dI ; dI/I

colour = colour_palette(0,4);

%% <dX^2>n vs In
Vart1 = Lum(:, 1:end-1);
xl = 'I_n';
Vart2 = dX_corr(:, 1:end).^2;
b = 18;
yl = '<\delta\theta^2_n>';

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, b, b+10);

s1 = std(v2bin,1,2);
b1 = binvals;

%***
f = figure;
errorbar(binvals, mean(v2bin, 2), std(v2bin,1,2)/sqrt(elts_per_bin), 'Color', colour(2,:), 'Linewidth', 1.5)
xlim([floor(min(binvals)*1.2) max(binvals)])
xticks([0:0.05:0.22])
xlabel(xl)
ylabel(yl)
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%% <dX^2>n &vs dIn-1
Vart1 = dLum(:, 1:end-1);
xl = '\DeltaI_n_-_1';
Vart2 = dX(:, 2:end).^2;
b = 17;
yl = '<\delta\theta^2>_n';

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, b, b+10);

s2 = std(v2bin,1,2);
b2 = binvals;

%***
f = figure;
errorbar(binvals,  mean(v2bin, 2), std(v2bin,1,2)/sqrt(elts_per_bin), 'k', 'Linewidth', 1.5)
xlim([-0.08 max(binvals)*1.2])
xticks([-0.08:0.04:0.08])
xlabel(xl)
ylabel(yl)
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';


%% <dX^2>n &vs dI/I(n-1)
colour = colour_palette(0,4);

Vart1_b = dLum(:, 1:mednbouts)./( (Lum(:, 1:mednbouts)+ Lum(:, 2:mednbouts+1))./2 ) ;
xl = '\deltaI/I_n_-_1';
Vart2 = dX(:, 2:mednbouts+1).^2;

b = 22;
yl = '<\delta\theta^2>_n';

[binvals, elts_per_bin, v2bin, ~, binedges] = BinsWithEqualNbofElements(Vart1_b, Vart2, b, b+3);
m = mean(v2bin,2);
s = std(v2bin, 1, 2);

%***
colour_spont = colour_palette(0,1);
f = figure;
plot([binvals(1)*1.2, binvals(end)], [0.18 0.18],...
    '--', 'Linewidth', 1.5, 'Color', colour_spont(1,:), 'DisplayName', 'spontaneous value')
hold on
errorbar(binvals,  m, s/sqrt(elts_per_bin), s/sqrt(elts_per_bin),...
    diff(binedges)/4, diff(binedges)/4,...
    'Linewidth', 1, 'Color', colour(2,:), 'Marker', 'sq', 'MarkerEdgeColor', colour(2,:), 'HandleVisibility', 'off')
xlim([min(binvals)*1.2 max(binvals)*1.2])
xticks([-1.5:0.5:1.5])
xlabel(xl)
ylabel(yl)
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
legend

% --- significance check ---
signif_val = 0.001;
[pm, pv] = mean_var_testing(Vart1_b,dX(:, 2:mednbouts+1), b, signif_val)

%% <dX^2>n &vs dI/I(n-1) per fish
Vart1_b = dLum(:, 1:end-1)./( (Lum(:, 1:end-2)+ Lum(:, 2:end-1))./2 ) ;
xl = '\deltaI/I_n_-_1';
Vart2 = dX_corr(:, 2:end).^2;

b = 7;
yl = '<\delta\theta^2>_n';

fish = unique(FishID);

beta = NaN(length(fish),1);
betaci = NaN(length(fish),2);
medboutsperfish = NaN(length(fish),1);
minboutsperfish = NaN(length(fish),1);
maxboutsperfish = NaN(length(fish),1);
meanXperfish = NaN(length(fish),1);
meanRonseqperfish = NaN(length(fish),1);
stdRonseqperfish = NaN(length(fish),1);
Rcircperfish = NaN(length(fish),1);

for i = 1:length(fish)
    seqafish = (FishID == i);
    var1 = Vart1_b(seqafish,:);
    var2 = Vart2(seqafish,:);
    if length(var1(:))-sum(isnan(var1(:))) < b
        continue
    end
    [binvals, elts_per_bin, v2bin, ~, binedges] = BinsWithEqualNbofElements(var1, var2, b, b+3);
    m = mean(v2bin,2);
    s = std(v2bin, 1, 2);
    medboutsperfish(i) = median(size(var1,2) - sum(isnan(var1),2));
    minboutsperfish(i) = min(size(var1,2) - sum(isnan(var1),2));
    maxboutsperfish(i) = max(size(var1,2) - sum(isnan(var1),2));
    linfittype = fittype('a*x','coefficients', {'a'});
    myfit = fit(binvals(binvals<0), m(binvals<0), linfittype);
    beta(i) = myfit.a;
    betaci(i,:) = confint(myfit)';
    plot(binvals, m);%, 'sq', 'MarkerFaceColor', 'k')
    hold on
    plot(binvals(binvals<0), myfit.a*binvals(binvals<0), 'LineWidth', 1.5)
    
    xfish = X(seqafish,2:mednbouts)+pi/2;
    meanRonseqperfish(i) = nanmean(circ_r(xfish'));
    stdRonseqperfish(i) = nanstd(circ_r(xfish'));
    xfish = xfish(:);
    xfish(isnan(xfish))=[];
    meanXperfish(i) = circ_mean(xfish);
    Rcircperfish(i) = circ_r(xfish);
end

figure;
errorbar(medboutsperfish, beta, beta-betaci(:,1), betaci(:,2)-beta, minboutsperfish, maxboutsperfish,...
    'o', 'LineWidth', 1.5, 'Color', 'k', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')


Rprojperfish = Rcircperfish.*cos(meanXperfish+pi/2);
meanRprojonseqperfish = meanRonseqperfish.*cos(meanXperfish+pi/2);

%***
figure
plot(beta, Rprojperfish, 'sq', 'MarkerSize', 8)
xlabel('standard deviation pente for negative \deltaI/I')
ylabel('<R.cos(\theta)>')
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Times New Roman';
 
%% <dX^2>n &vs dI/I(n-1) with dt
b = 12;
f = figure;
plot([binvals(1), binvals(end)], [0.18 0.18],...
    '--', 'Linewidth', 1.5, 'Color', colour_spont(1,:), 'DisplayName', 'spontaneous')
hold on
for dt = 0 : 3
    dII = dLum(:, 1:end-(1+dt))./( (Lum(:, 1:end-(2+dt)) + Lum(:, 2:end-(1+dt)))./2 ) ;
    dXsq = dX(:, 2+dt:end).^2;
    [binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(dII, dXsq, b, b+10);
    
    %***
    errorbar(binvals,  mean(v2bin, 2), std(v2bin,1,2)/sqrt(elts_per_bin),...
        'Linewidth', 1.2, 'Color', colour(dt+2,:), 'DisplayName', ['dn = ' num2str(dt+1)])
    hold on
end

xl = 'dI/I_{n-dn}';
yl = '<\delta\theta^2_n>';
xlabel(xl)
ylabel(yl)
legend
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
