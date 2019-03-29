function[f] = half2gauss(x, a, muf, mut, sigmaf, sigmat)

muf = 0;
sigmaf = 0.1;

f = (1-a)/(sigmaf*sqrt(2*pi))* exp(-(x-muf).^2/(2*sigmaf^2)) + a/(sigmat*sqrt(2*pi))* exp(-(x-mut).^2/(2*sigmat^2));

end

