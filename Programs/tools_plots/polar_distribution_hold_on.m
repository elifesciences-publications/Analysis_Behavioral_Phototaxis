function[fig] = polar_distribution_hold_on(x, Nbins, FishID, selected_color)


%% appearance
[colour] = selected_color;

%%
% --- global distribution ---
angle = ((1:Nbins) - 0.5) * 2*pi/Nbins;
xpdf = hist(mod(x(:),2*pi),Nbins)/sum(hist(mod(x(:),2*pi),Nbins));

circr = circ_r(wrapToPi(x(:)));
theta_mean = circ_mean(wrapToPi(x(:)));
Rproj_mean = circ_r(x(:)).*cos(theta_mean);

% --- stats on individual fish ---
different_fish = unique(FishID);
xpdf_fish = NaN(different_fish(end), Nbins);
for i = 1 : different_fish(end)
    seqs = i==FishID;
    xf = x(seqs,:);
    xpdf_fish(i,:) = hist(mod(xf(:),2*pi),Nbins)/sum(hist(mod(xf(:),2*pi),Nbins));
end
xmean_fish = nanmean(xpdf_fish,1);
xstd_fish = nanstd(xpdf_fish,1,1);

%***

% --- fake plot to have the same scale ---
fake = polar(0, 0.15);
fake.HandleVisibility = 'off';
hold on

% --- plot the distribution ---
polarwitherrorbar([angle, angle(1)], smooth([xmean_fish, xmean_fish(1)])',...
    [ xstd_fish,  xstd_fish(1)]/sqrt(different_fish(end)), colour);
hold on
p1 = polar([angle, angle(1)],smooth([xpdf, xpdf(1)])');
p1.LineStyle = '--';
p1.LineWidth = 2;
p1.Color =  colour;
p1.DisplayName = ['<\theta> = ' num2str(theta_mean, 1)];

% --- plot and display the circular mean and R ---
hold on
p2=polar([0 theta_mean], [0 circr/10]);
p2.DisplayName = ['<R.cos(\theta)> = ' num2str(Rproj_mean, 2)];
p2.LineWidth = 2;
p2.Color =  colour;

% --- polar plot parameters ---
theta_ticks_to_remove = {'30','60','120','150','210','240','300','330'};
for th = 1:length(theta_ticks_to_remove)
    set(findall(gcf, 'String', theta_ticks_to_remove{th}) ,'String', ' ');
end
theta_ticks = {'90', '180', '270'};
theta_labels = {'-\pi/2', '\pi', '\pi/2'};
for th = 1:length(theta_ticks)
    set(findall(gcf, 'String', theta_ticks{th}),...
        'String', theta_labels{th}, 'FontSize', 16, 'FontName', 'Times New Roman');
end

% Altering the radial label
%set(findall(gcf, 'String', '  4'),'String', ' Four'); 

ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'normal';
view([-90 -90])