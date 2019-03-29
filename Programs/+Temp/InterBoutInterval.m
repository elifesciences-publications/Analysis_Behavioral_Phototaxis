%% Inter-bout interval in temporal phototaxis data

%% IBI vs I & dI/I
Vart1_a = Lum(:, 1:end-2);
Vart1_b = dLum(:, 1:end-1)./( (Lum(:, 1:end-2) + Lum(:, 2:end-1))./2 ) ;
Vart2 = IBI(:, 2:end);

xla = 'I';
xlb = 'dI/I';
b = 14;
yl = 'bout frequency (Hz)';

[binvals_a, elts_per_bin, v2binMatrix_a] = BinsWithEqualNbofElements(Vart1_a, Vart2, b, b+10);
[binvals_b, elts_per_bin, v2binMatrix_b] = BinsWithEqualNbofElements(Vart1_b, Vart2, b, b+10);

%***
f = figure;
subplot(1,2,1)
errorbar(binvals_a,  1./mean(v2binMatrix_a, 2), 1./std(v2binMatrix_a,1,2)/sqrt(elts_per_bin),...
    'k', 'LineWidth', 1.5)
xlabel(xla)
ylabel(yl)
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

subplot(1,2,2)
errorbar(binvals_b,  1./mean(v2binMatrix_b, 2), 1./std(v2binMatrix_b,1,2)/sqrt(elts_per_bin),...
    'k', 'LineWidth', 1.5)
xlabel(xlb)
ylabel(yl)
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';