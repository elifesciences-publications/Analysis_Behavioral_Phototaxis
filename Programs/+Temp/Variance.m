%% Variance analysis of temporal phototaxis data
% dXsq versus I ; dI ; dI/I


%% <dX^2>n vs In
Vart1 = Lum(:, 1:end-1);
xl = 'I_n';
Vart2 = dX_corr(:, 1:end).^2;
b = 18;
yl = '<\delta\theta^2_n>';

[binvals, elts_per_bin, v2bin_pos] = BinsWithEqualNbofElements(Vart1, Vart2, b, b+10);

s1 = std(v2bin_pos,1,2);
b1 = binvals;

%***
f = figure;
errorbar(binvals, mean(v2bin_pos, 2), std(v2bin_pos,1,2)/sqrt(elts_per_bin), 'Color', colour(2,:), 'Linewidth', 1.5)
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

[binvals, elts_per_bin, v2bin_pos] = BinsWithEqualNbofElements(Vart1, Vart2, b, b+10);

s2 = std(v2bin_pos,1,2);
b2 = binvals;

%***
f = figure;
errorbar(binvals,  mean(v2bin_pos, 2), std(v2bin_pos,1,2)/sqrt(elts_per_bin), 'k', 'Linewidth', 1.5)
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
Vart2(Vart2>50) = NaN;
b = 22;
yl = '<\delta\theta^2>_n';

[binvals, elts_per_bin, v2bin_pos, ~, binedges] = BinsWithEqualNbofElements(Vart1_b, Vart2, b, b+10);

s3 = std(v2bin_pos,1,2);
b3 = binvals;

%***
colour_spont = colour_palette(0,1);
f = figure;
plot([binvals(1)*1.2, binvals(end)], [0.18 0.18],...
    '--', 'Linewidth', 1.5, 'Color', colour_spont(1,:), 'DisplayName', 'spontaneous value')
hold on
errorbar(binvals,  mean(v2bin_pos, 2), std(v2bin_pos,1,2)/sqrt(elts_per_bin), std(v2bin_pos,1,2)/sqrt(elts_per_bin),...
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

%% <dX^2>n &vs dI/I(n-1) with dt
b = 12;
f = figure;
plot([binvals(1), binvals(end)], [0.18 0.18],...
    '--', 'Linewidth', 1.5, 'Color', colour_spont(1,:), 'DisplayName', 'spontaneous')
hold on
for dt = 0 : 3
    dII = dLum(:, 1:end-(1+dt))./( (Lum(:, 1:end-(2+dt)) + Lum(:, 2:end-(1+dt)))./2 ) ;
    dXsq = dX(:, 2+dt:end).^2;
    [binvals, elts_per_bin, v2bin_pos] = BinsWithEqualNbofElements(dII, dXsq, b, b+10);
    
    %***
    errorbar(binvals,  mean(v2bin_pos, 2), std(v2bin_pos,1,2)/sqrt(elts_per_bin),...
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
