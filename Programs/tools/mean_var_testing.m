function[pm, pv] = mean_var_testing(v1, v2, b, signif_val)

% Variance and mean analysis tests

Var1 = v1 ;
Var2 = v2;
Var2sq = v2.^2;

[binvals1, elts_per_bin1, v2bin1, ~, binedges1] = BinsWithEqualNbofElements(Var1, Var2, b, b+3);
m = mean(v2bin1, 2);
[binvals2, elts_per_bin2, v2bin2, ~, binedges2] = BinsWithEqualNbofElements(Var1, Var2sq, b, b+3);
s = mean(v2bin2, 2);

% --- transform to perform mean and variance tests ---
vtest = v2bin1';

% - test on mean -
% ttest H(0) : data tested comes from a normal distribution with mean = 0 and unknown variance
[hm, pm] = ttest(vtest,0, signif_val)
% for plot
xmeantest = binvals1;
ymeantest = 1.2*hm.*m';
xmeantest(hm==0) = [];
ymeantest(hm==0) = [];

% - test of variance -
bin0 = find(binvals1==0);
if isempty(bin0)
    [~, bin0] = min(abs(binvals1));
end
v2bins0 = v2bin1(bin0,:);
ref_var = var(v2bins0(:));
[hv, pv] = vartest(vtest, ref_var, 'Alpha', signif_val)
% for plot
xvartest = binvals2;
xvartest(hv==0) = [];
yvartest = 1.2*hv.*s';
yvartest(hv==0) = [];

% - 2-by-2 variance test -
% - test of variance -
vtest1 = vtest(:,1:end-1);
vtest2 = vtest(:,2:end);
[hv2, pv2] = vartest2(vtest1, vtest2, 'Alpha', signif_val)
% for plot
xvartest2 = binvals2(1:end-1) + diff(binvals2)/2;
xvartest2(hv2==0) = [];
yvartest2 = 1.2*hv2.*s(1:end-1)';
yvartest2(hv2==0) = [];

%***
figure
errorbar(binvals1,  mean(v2bin1, 2), std(v2bin1,1,2)/sqrt(elts_per_bin1), std(v2bin1,1,2)/sqrt(elts_per_bin1),...
    diff(binedges1)/4, diff(binedges1)/4)
hold on
plot(xmeantest, ymeantest, '*')
title(['mean significantly =/= 0. Significance level : ' num2str(signif_val)])


figure
errorbar(binvals2,  mean(v2bin2, 2), std(v2bin2,1,2)/sqrt(elts_per_bin2), std(v2bin2,1,2)/sqrt(elts_per_bin2),...
    diff(binedges2)/4, diff(binedges2)/4)
hold on
plot(xvartest, yvartest, '*')
title(['var significantly =/= reference value ' num2str(ref_var) '. Significance level : ' num2str(signif_val)])

figure
errorbar(binvals2,  mean(v2bin2, 2), std(v2bin2,1,2)/sqrt(elts_per_bin2), std(v2bin2,1,2)/sqrt(elts_per_bin2),...
    diff(binedges2)/4, diff(binedges2)/4)
hold on
plot(xvartest2, yvartest2, '*')
title(['consecutive vars significantly =/= from one another. Significance level : ' num2str(signif_val)])
