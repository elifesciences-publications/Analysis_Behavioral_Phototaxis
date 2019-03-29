function[wturn, wfor, pturn, pfor, cipturn] = Gauss2custom(v2binMatrix)

elts_per_bin = size(v2binMatrix,2);
minv2 = floor(prctile(v2binMatrix(:), 0.1)*10)/10;
maxv2 = ceil(prctile(v2binMatrix(:), 99.9)*10)/10;

% custom fit on bins
mabsdX = nanmean(abs(v2binMatrix),2);
wsq = nanvar(v2binMatrix, 1, 2);

Nbins = size(v2binMatrix,1);
wturn = NaN(1, Nbins);
wfor = NaN(1, Nbins);
pfor = NaN(1, Nbins);
pturn = NaN(1, Nbins);
cipturn = NaN(Nbins, 2);

for i = 1 : Nbins
    disp(i)
    
    binwidth = 4*iqr(v2binMatrix(i,:))/((elts_per_bin)^(1/3));
    bins = minv2 : binwidth : maxv2;
    [~, devfromzeroidx] = min(abs(bins));
    bins = bins -  bins(devfromzeroidx);
    [d_hist, x_hist]  = histcounts(v2binMatrix(i,:), bins);
    x_hist = x_hist(1:end-1)+ mean(binwidth);
    d_hist = d_hist/(sum(d_hist)*binwidth);
    
    %constrained fit (imposed variance and mean of the abs)
    alpha = mabsdX(i);
    beta = wsq(i);
    
    ft = fittype('double_gauss_constrained(x,a1,alpha,beta)', 'problem', {'alpha','beta'});
    f=fit(x_hist', d_hist', ft, 'StartPoint', [0.4], 'Lower', 0, 'Upper', 1, 'problem', {alpha,beta});
    a = f.a1;
    wturn(i) = (alpha*a + sqrt(alpha^2*a^2 - (alpha^2 - beta*(1-a))*(a^2+a*(1-a))))/(a^2+a*(1-a));
    wfor(i) = (alpha*sqrt(pi/2) - a*wturn(i))/(1-a);
    pfor(i) = 1-a;
    pturn(i) = a;
    cipturn(i,:) = confint(f)';
     
    h=figure;
    hold off; plot(x_hist, d_hist);
    hold on;  plot(f, x_hist',d_hist')
    waitfor(h)  
end

