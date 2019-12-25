function[dt2] = amplitude_fit(x, wturn, wfor, aturn, k2)

wt=wturn^2;
wf=wfor^2;
 
afor=1-aturn;

f = afor*normpdf(x,0,wfor);
g = aturn*normpdf(x,0,wturn);
h = g./(f+g);

alpha=aturn/(1-aturn);
vardt=aturn*wt+(1-aturn)*wf;

dt2=wf+k2*alpha*(wt-wf)+h*(wt-wf)*(1-k2-k2*alpha);

dt2 = dt2./vardt;

end