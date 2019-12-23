%% Autocorrelations in temporal phototaxis data

%% 1-Q
Vart1 = dLum(:, 1:end-2)./( (Lum(:, 1:end-3) + Lum(:, 2:end-2))./2 ) ;
xl = 'dI/I';
Vart2 = dX(:, 2:end-1).*dX(:, 3:end)./(abs(dX(:, 2:end-1)).*abs(dX(:, 3:end)));
b = 10;
yl = 'autocorr';

[binvals, elts_per_bin, v2bin_pos] = BinsWithEqualNbofElements(Vart1, Vart2, b, b+10);

%***
f = figure;
errorbar(binvals,  mean(v2bin_pos, 2), std(v2bin_pos,1,2)/sqrt(elts_per_bin), 'k')
xlabel(xl)
ylabel(yl)


%% spot autocorrelations

[dimdLum_idx1, dimdLum_idx2] = find(dLum(:, 1:end-1) < 0);
[augdLum_idx1, augdLum_idx2] = find(dLum(:, 1:end-1) > 0);

VartA1 = NaN(size(dX));
VartA2 = NaN(size(dX));
VartB1 = NaN(size(dX));
VartB2 = NaN(size(dX));
VartA1(dimdLum_idx1, dimdLum_idx2) = dX(dimdLum_idx1, dimdLum_idx2);
VartA2(dimdLum_idx1, dimdLum_idx2) = dX(dimdLum_idx1, dimdLum_idx2+1);

VartB1(augdLum_idx1, augdLum_idx2) = dX(augdLum_idx1, augdLum_idx2);
VartB2(augdLum_idx1, augdLum_idx2) = dX(augdLum_idx1, augdLum_idx2+1);

[binvalsA, elts_per_binA, v2AbinMatrix] = BinsWithEqualNbofElements(VartA1, VartA2, 10, 20);
[binvalsB, elts_per_binB, v2BbinMatrix] = BinsWithEqualNbofElements(VartB1, VartB2, 10, 20);

errorbar(binvalsA,  mean(v2AbinMatrix, 2), std(v2AbinMatrix,1,2)/sqrt(elts_per_binA), 'b')
hold on
errorbar(binvalsB,  mean(v2BbinMatrix, 2), std(v2BbinMatrix,1,2)/sqrt(elts_per_binB), 'k')
legend('neg', 'pos')

%% dX(n)  & dX(n-1)
dII_n = dLum(:, 1:end-2)./( (Lum(:, 1:end-3) + Lum(:, 2:end-2))./2 ) ;
dXn =  dX_corr(:, 1:end-2);
dXnp1 =  dX_corr(:, 2:end-1);
dXndXnp1 = dXn.*dXnp1./(abs(dXn).*abs(dXnp1));
b=10;

thresh = 0.1;
posddI = dII_n > thresh;
negddI = dII_n <-thresh;

[binvals_pos, elts_per_bin_pos, v2bin_pos, ~] = BinsWithEqualNbofElements(dXn(posddI), dXnp1(posddI), b, b+3);
[binvals_neg, elts_per_bin_neg, v2bin_neg, ~] = BinsWithEqualNbofElements(dXn(negddI), dXnp1(negddI), b, b+3);

%***
f = figure;
hold on
errorbar(binvals_pos,  mean(v2bin_pos, 2), std(v2bin_pos,1,2)/sqrt(elts_per_bin_pos),...
    'Color', colour(3,:), 'LineWidth', 1.5, 'DisplayName', 'dII_n > threshold')
errorbar(binvals_neg,  mean(v2bin_neg, 2), std(v2bin_neg,1,2)/sqrt(elts_per_bin_neg),...
    'Color', colour(4,:), 'LineWidth', 1.5, 'DisplayName', 'dII_n < -threshold')
legend
xlim([-pi/2 pi/2])
xticks([-pi/2,-pi/4, 0,pi/4, pi/2])
xticklabels({'-\pi/2','-\pi/4', '0','\pi/4', '\pi/2'})
ylim([-0.4 0.4])
ax=gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 14;

%%
dII = dLum(:, 1:end-1)./( (Lum(:, 1:end-2)+ Lum(:, 2:end-1))./2 ) ;
bins = 12;
<<<<<<< HEAD
Temp.auto_co_reinforcement(dX, dII, bins);
=======
Temp.a.auto_co_reinforcement(dX, dII, bins);
>>>>>>> master
