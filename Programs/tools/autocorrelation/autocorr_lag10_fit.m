function f = autocorr_lag10_fit(x, p_flip, aturn, wturn, wfor)

p_TF = 1-aturn;

f=2/pi*(1-2*p_flip).^x*((1-p_TF)*wturn)^2/((1-p_TF)*wturn^2+p_TF*wfor^2);
f(1)=1;

end