function[p] = two_sample_distribution_test(x, xs)

med_bouts_per_seq = median(size(x,2)-sum(isnan(x),2));


x_samp = x(:,2:med_bouts_per_seq);
x_samp = x_samp(:);
nans_ind = isnan(x_samp);
x_samp(nans_ind)=[];

xs_samp = xs(1:size(x,1),2:med_bouts_per_seq);
xs_samp = xs_samp(:);
xs_samp(nans_ind) = [];

x1 = wrapToPi(x_samp);
x2 = wrapToPi(xs_samp);

histogram(x1)
hold on
histogram(x2)

[p,U2]=circ_kuipertest(x1,x2)
[pval med P] = circ_cmtest(x1, x2)

% non...
%[p, k, K] = circ_kuipertest(x_samp, xs_samp, 1000, 1)

end