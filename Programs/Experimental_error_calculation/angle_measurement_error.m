function [a, mu, sigma, se] = angle_measurement_error(angle_for_error, boutIndices_ini, framerate, fig)

tau = ceil(0.2*framerate);

Baseline = NaN(numel(angle_for_error),1);
count = 1;
for i = 1 : size(angle_for_error,1)
    afe = angle_for_error(i,:);
    baseline = afe' - smooth(afe, framerate(i)/2);
    bi = boutIndices_ini(i,:);
    for j = 1 : length(bi)-sum(isnan(bi))
        if bi(j)>tau(i) 
            baseline(bi(j)-tau(i) : bi(j)+tau(i)) = nan;
        else
            baseline(1 : bi(j)+tau(i)) = nan;
        end        
    end
    baseline(isnan(baseline)) = [];
    Baseline(count:count+(length(baseline))-1) = baseline;
    count = find(isnan(Baseline),1);
end

Baseline = deg2rad(Baseline);
se = nanstd(Baseline);

[y, x] = histcounts(Baseline,round(sqrt(length(Baseline))));
% --- matlab built-in gaussian fit ---
%f = fit(x',y','gauss1');
% if fig
%     plot(f, x, y);
%     hold on
%     xlim([-0.3 0.3])
% end
% a = f.a1;
% mu = f.b1;
% sigma = f.c1;

% --- another gaussian fit ---
binwidth = mean(diff(x));
x = x(1:end-1)+ binwidth;
y = y/(sum(y)*binwidth);
[sigma, mu] = gaussfit( x, y);
a = 1/(sqrt(2*pi)*sigma);

if fig
    plot(x,y)
    hold on
    plot(x, a*exp( -(x - mu).^2 / (2*sigma^2)))
    drawnow
end

