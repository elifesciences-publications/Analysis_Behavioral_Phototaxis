function[f, xhist, yhist, muabs, w, Wturn, Wfor, pturn] = Gauss2custom_simple(data)

data = data(:);

%constrained fit (imposed variance and mean of the abs)
muabs=nanmean(abs(data));
w=nanvar(data);

binwidth = 3*iqr(data)/((length(data)-sum(isnan(data)))^(1/3));
bins = min(data):binwidth:max(data);
[yhist, edges]  = histcounts(data, bins);
xhist=edges(1:end-1)+binwidth/2;
yhist=yhist/(sum(yhist)*binwidth);

ft = fittype('double_gauss_constrained(x, a1, muabs, w)', 'problem', {'muabs','w'});
f = fit(xhist', yhist', ft, 'StartPoint', [0.2], 'Lower', 0, 'Upper', 1, 'problem', {muabs,w}, 'Robust', 'LAR');
a = f.a1;
aci = confint(f,0.99);

wturn = (muabs*a+sqrt(muabs^2*a^2-(muabs^2-w*(1-a))*(a^2+a*(1-a))))/(a^2+a*(1-a));
wturnci1 = (muabs*aci(1)+sqrt(muabs^2*aci(1)^2-(muabs^2-w*(1-aci(1)))*(aci(1)^2+aci(1)*(1-aci(1)))))...
    /(aci(1)^2+aci(1)*(1-aci(1)));
wturnci2 = (muabs*aci(2)+sqrt(muabs^2*aci(2)^2-(muabs^2-w*(1-aci(2)))*(aci(2)^2+aci(2)*(1-aci(2)))))...
    /(aci(2)^2+aci(2)*(1-aci(2)));
Wturn = [wturn ; wturnci1 ; wturnci2];

wfor = (muabs*sqrt(pi/2)-a*wturn)/(1-a);
wforci1 = (muabs*sqrt(pi/2)-aci(1)*wturn)/(1-aci(1));
wforci2 = (muabs*sqrt(pi/2)-aci(2)*wturn)/(1-aci(2));
Wfor = [wfor ; wforci1 ; wforci2];

pturn = [a ; aci(1) ; aci(2) ]; 

%pfor = [1-a; 1-aci(2) ; 1-aci(1) ];



