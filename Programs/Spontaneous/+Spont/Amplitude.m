function[fig] = Amplitude(dXi)

% Check amplitude correlation
adXi = (dXi).^2;
adXn = adXi(:,1:end-1);
adXnp1 = adXi(:,2:end);

var1 = adXn;
var2 = adXnp1/nanvar(adXn(:));
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, 18, 24);
mv2 = mean(v2bin,2);
stdv2 = std(v2bin,1,2);

%*** 
fig = figure;

%set color palette
[colour] = colour_palette(0,1);

errorbar(binvals, mv2, stdv2/2,...
    'Color', colour(4,:), 'LineWidth', 1, 'DisplayName', '\pm \sigma/2')
hold on
errorbar(binvals, mv2, stdv2/sqrt(elts_per_bin-1),...
    'Color', colour(2,:), 'LineWidth', 1.5, 'DisplayName', '\pm \sigma/\surd{n-1}')

legend

xlabel('<|d\theta|>_n')
ylabel('<|d\theta|>_{n+1}')

xlim([0 1.5])
ylim([0 1])