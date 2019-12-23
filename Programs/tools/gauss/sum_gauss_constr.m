function[f] = sum_gauss_constr(x, a, mut, sigmat, muf, sigmaf)

f = (1-a)/(sigmaf*sqrt(2*pi))* exp(-(x-muf).^2/(2*sigmaf^2)) + a/(sigmat*sqrt(2*pi))* exp(-(x-mut).^2/(2*sigmat^2));

end

