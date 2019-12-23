function[fig] = polar_distribution_simu_exp(thetaData, thetaSimu, Nbins, FishID, colourID)


%% appearance
[colour] = colour_palette(0, colourID);

%%
% --- global distribution ---
angle = ((1:Nbins) - 0.5) * 2*pi/Nbins;
binsize = mean(diff(angle));
Nbinssimu = Nbins*8;
angle_simu = ((1:Nbinssimu) - 0.5) * 2*pi/Nbinssimu;
binsize_simu = mean(diff(angle_simu));
thetapdf = hist(mod(thetaSimu(:),2*pi),Nbinssimu)/(sum(hist(mod(thetaSimu(:),2*pi),Nbinssimu))*binsize_simu);

circr = circ_r(wrapToPi(thetaSimu(:)));
theta_mean_simu = circ_mean(wrapToPi(thetaSimu(:)));
Rproj_mean_simu = circ_r(thetaSimu(:)).*cos(theta_mean_simu);

theta_mean_data = circ_mean(wrapToPi(thetaData(:)));
Rproj_mean_data = circ_r(thetaData(:)).*cos(theta_mean_data);

% --- stats on individual fish ---
different_fish = unique(FishID);
thetapdf_fish = NaN(different_fish(end), Nbins);
for i = 1 : different_fish(end)
    seqs = i==FishID;
    thetaf = thetaData(seqs,:);
    thetapdf_fish(i,:) = hist(mod(thetaf(:),2*pi),Nbins)/(sum(hist(mod(thetaf(:),2*pi),Nbins))*binsize);
%     polarplot([angle, angle(1)],[xpdf_fish(i,:), xpdf_fish(i,1)])
%     hold on
end
theta_mean_fish = nanmean(thetapdf_fish,1);
theta_std_fish = nanstd(thetapdf_fish,1,1);

%***
fig = figure;

% --- fake plot to have the same scale ---
fake = polar(0, 0.3);
fake.HandleVisibility = 'off';
hold on

% --- plot the distribution ---
polarwitherrorbar([angle, angle(1)], smooth([theta_mean_fish, theta_mean_fish(1)])',...
    [ theta_std_fish,  theta_std_fish(1)]/sqrt(different_fish(end)), colour(5,:));
hold on
p1 = polar([angle_simu, angle_simu(1)], [thetapdf, thetapdf(1)]);
p1.LineWidth = 2;
p1.Color =  colour(2,:);
p1.DisplayName = ['<\theta>_{simu} = ' num2str(wrapToPi(theta_mean_simu),2) ' rad'...
    '(data : ' num2str(wrapToPi(theta_mean_data),2) ')'];

% --- plot and display the circular mean and R ---
hold on
p2=polar([0 theta_mean_simu], [0 circr/10]);
p2.DisplayName = ['<R_{simu}> = ' num2str(Rproj_mean_simu,3)...
    '(data : ' num2str(Rproj_mean_data,3) ')'];
p2.LineWidth = 2;
p2.Color =  colour(3,:);

% --- polar plot parameters ---
theta_ticks_to_remove = {'30','60','120','150','210','240','300','330'};
for th = 1:length(theta_ticks_to_remove)
    set(findall(gcf, 'String', theta_ticks_to_remove{th}) ,'String', ' ');
end
theta_ticks = {'90', '180', '270'};
theta_labels = {'-\pi/2', '\pi', '\pi/2'};
for th = 1:length(theta_ticks)
    set(findall(gcf, 'String', theta_ticks{th}),...
        'String', theta_labels{th}, 'FontSize', 14, 'FontName', 'Times New Roman');
end

% Altering the radial label
%set(findall(gcf, 'String', '  4'),'String', ' Four'); 

ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'normal';
view([-90 -90])

