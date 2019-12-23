%% error  fit extr

[~, mu, sigma] = theta_measure_error_estimation(0);
%%
[colour] = colour_palette(0,1);

% dX vs X
minbin = 9;
maxbin = 12;

Vart1 = Xiwrapped(:,1:end-1);
Vart2 = dXi;
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, minbin, maxbin);

%***
fig = figure;
hold on

%***
area([-pi pi], [sigma/2 sigma/2], -sigma/2,...
    'FaceColor', colour(5,:), 'FaceAlpha', 0.3, 'DisplayName', 'detection error')

%***
plot([min(Xiwrapped(:)), max(Xiwrapped(:))], [0 0],...
    '--', 'Color', [0.7 0.7 0.7], 'DisplayName', 'y=0');
errorbar(binvals, mean(v2bin,2), std(v2bin,1,2)/sqrt(elts_per_bin-1),...
    'LineWidth', 1.5, 'Color', colour(5,:), 'DisplayName', '<\delta\theta>')

% -----
Vart2 = dXissbiais;
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, minbin,maxbin);

%***
errorbar(binvals, mean(v2bin,2), std(v2bin,1,2)/sqrt(elts_per_bin-1),...
    'LineWidth', 1.5, 'Color', colour(3,:), 'DisplayName',  '<\delta\theta> - <\delta\theta>_f_i_s_h')

xlabel('\theta (rad)')
ylabel('<\delta\theta> (rad)')
ylim([-0.1 0.1])
xlim([-pi pi])
xticks([-pi, -pi/2, 0, pi/2, pi])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
legend
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
