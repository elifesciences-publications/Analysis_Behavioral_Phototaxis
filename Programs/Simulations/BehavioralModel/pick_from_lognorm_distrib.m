function[projection] = pick_from_lognorm_distrib(x, sigma, mu, N)

% input
% possible values of x = 0:0.01:30;
% sigma = 0.6;
% mu=0.8;


lognorm_distribution = @(x) 1./(x*sigma*sqrt(2*pi)) .* exp(-(log(x)-mu).^2/(2*sigma^2));
cdf = cumsum(lognorm_distribution(x), 'omitnan')/nansum(lognorm_distribution(x));

[cdf, mask] = unique(cdf);
x = x(mask);

randomValues = rand(1, N);

% inverse interpolation to achieve P(x) -> x projection of the random values
projection = interp1(cdf, x, randomValues);