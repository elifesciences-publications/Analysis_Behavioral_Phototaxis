function testBiasedDiffModel ()


x = 0:0.01:pi;
a = mean(A);
d = 0.26;

k = a/d;

N =  2*k/(1-exp(-2*k*pi));
ps = N*exp(-2*k*x);

figure
h = histogram(Xwrapped, 'Normalization', 'count');
maxcount = max(h.Values);
hold on
plot(-x+pi, maxcount*ps/max(ps))
ylim([0 maxcount])
