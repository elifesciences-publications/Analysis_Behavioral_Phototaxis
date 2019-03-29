function f = autocorr_lag1_fit(x, p_flip, pturn, wturn, wfor)

pfor=1-pturn;   % ratio of forward bouts
p_TF = 1-pturn;   % probability of going from turn to forward state

f = pfor*normpdf(x,0,wfor);
g = pturn*normpdf(x,0,wturn);
h = g./(f+g);

%f = sqrt(2/pi)*(1-2*p_flip)*sign(x).*(h*(1-p_TF)*wturn);
f = (1-2*p_flip)*sign(x).*(h*(1-p_TF)*wturn);

end