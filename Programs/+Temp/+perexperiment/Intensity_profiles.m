
%% exponential complete

[colour] = temp_colours();


experiment = cell(3,1);
experiment{1} = 'exp60';
experiment{2} = 'exp30';
experiment{3} = 'sin60';

experiment_txt{1} = 'Exponential 60%';
experiment_txt{2} = 'Exponential 30%';
experiment_txt{3} = 'Sinusoidal 60%';

savefig = 1;
for i = 1 : 3
    [~, ~, ~, lum_th, ~] = Temp.chooseExpType(experiment{i});
    
    % luminosity
    theta_rad = deg2rad(lum_th(:,1));
    lumW = lum_th(:,2);
    
    %***
    fig = figure;
    fig.Name = 'Intensity profile';
    
    s = subplot(2,2,1);
    plot(theta_rad, lumW,'k', 'LineWidth', 2);
    xlabel('fish orientation (rad)');
    ylabel('I (W/cm^2)');
    xticks([0 pi/4 pi/2 3*pi/4 pi])
    xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
    xlim([0 pi])
    title('')
    set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
    
    subplot(2,2,2)
    plot(lumW(2:end), diff(lumW),'k', 'LineWidth', 2);
    xlabel('I (W/cm^2)');
    ylabel('\Delta I');
    xlim([0 lumW(end)])
    set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
    
    subplot(2,2,3:4)
    plot(theta_rad(2:end), diff(lumW)./lumW(2:end),'k', 'LineWidth', 2)
    xticks([0 pi/4 pi/2 3*pi/4 pi])
    xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
    xlabel('orientation (rad)');
    ylabel('\DeltaI/I');
    xlim([0 pi])
    set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
    
    title(s, experiment_txt{i})
end

%% exponential : only basic profiles
[colour] = Temp.temp_colours();

experiment = cell(3,1);
experiment{1} = 'sin60';
experiment{2} = 'exp60';
experiment{3} = 'exp30';

experiment_txt{1} = 'Sinusoidal 60%';
experiment_txt{2} = 'Exponential 60%';
experiment_txt{3} = 'Exponential 30%';
%***
fig = figure;
fig.Name = 'Intensity profile';
hold on

for i = 1 : 3
    [~, ~, ~, lum_th, ~] = Temp.chooseExpType(experiment{i});
    
    % luminosity
    theta_rad = abs(deg2rad(lum_th(:,1))-pi);
    lumW = lum_th(:,2);    
    plot(theta_rad, lumW, 'DisplayName', experiment{i},...
        'LineWidth', 2, 'Color', colour(i,:));
    xlabel('fish orientation (rad)');
    ylabel('I (\muW/cm^2)');
end

xticks([0 pi/2 pi])
xticklabels({'0','\pi/2','\pi'})
xlim([0 pi])
ylim([0 300])
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
title('Intensity profiles')
legend

%% exponential derivative all on one figure

experiment = cell(3,1);
experiment{1} = 'sin60';
experiment{2} = 'exp60';
experiment{3} = 'exp30';

experiment_txt{1} = 'Sinusoidal 60%';
experiment_txt{2} = 'Exponential 60%';
experiment_txt{3} = 'Exponential 30%';

savefig = 0;
%***
fig = figure;
fig.Name = 'All temporal intensity profiles';
colors = [ 145/255 171/255 60/255; 229/255 98/255 38/255; 89/255 10/255 49/255];
hold on
for i = 1 : 3
    [~, ~, ~, lum_th, ~] = Temp.chooseExpType(experiment{i});
    color = colors(i,:);
    
    % luminosity
    theta_rad = deg2rad(lum_th(:,1));
    lumW = lum_th(:,2);
    
 %   subplot(2,2,1);
%     hold on
%     plot(theta_rad, lumW, 'LineWidth', 2, 'Color', color);
%     xlabel('fish orientation (rad)');
%     ylabel('I (W/cm^2)');
%     xticks([0 pi/4 pi/2 3*pi/4 pi])
%     xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
%     xlim([0 pi])
%     title('')
%     set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
    
    %subplot(2,2,2)
    subplot(2,1,1)
    hold on
    plot(lumW(2:end), diff(lumW), 'LineWidth', 2, 'Color', color);
    xlabel('I (W/cm^2)');
    ylabel('\Delta I');
    xlim([0 lumW(end)])
    set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
    
    %subplot(2,2,3:4)
    subplot(2,1,2)
    hold on
    plot(theta_rad(2:end), log(diff(lumW)./lumW(2:end)), 'LineWidth', 2, 'Color', color)
    xticks([0 pi/4 pi/2 3*pi/4 pi])
    xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
    xlabel('orientation (rad)');
    ylabel('log(\DeltaI/I)');
    xlim([0 pi])
    set(gca, 'FontSize', 14, 'FontName', 'Times New Roman')
end
legend(experiment_txt)

