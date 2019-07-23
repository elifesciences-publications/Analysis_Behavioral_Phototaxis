function[fig] = polar_distribution(x, Nbins, FishID, colourID)


%% appearance
[colour] = colour_palette(0, colourID);
txtAngle = [0,0];
txtRng = [-0.08, -0.06];

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
%     polarplot([angle, angle(1)],[xpdf_fish(i,:), xpdf_fish(i,1)])
%     hold on
end
xmean_fish = nanmean(xpdf_fish,1);
xstd_fish = nanstd(xpdf_fish,1,1);

% --- gen text for plot ---
names = {['<\theta> = ' num2str(round(wrapToPi(theta_mean)*100)/100) ' rad'],...
    ['<R> = ' num2str(round(Rproj_mean*1000)/1000)]};

%***
fig = figure;

% --- fake plot to have the same scale ---
fake = polar(0, 0.25);
fake.HandleVisibility = 'off';
hold on

% --- plot the distribution ---
polarwitherrorbar([angle, angle(1)], smooth([xmean_fish, xmean_fish(1)])',...
    [ xstd_fish,  xstd_fish(1)]/sqrt(different_fish(end)), colour(5,:));
hold on
%p1 = polar([angle, angle(1)],smooth([xpdf, xpdf(1)])');
%p1.LineWidth = 2;
%p1.Color =  colour(2,:);
%p1.DisplayName = '<\theta>';

% --- plot and display the circular mean and R ---
hold on
p2=polar([0 theta_mean], [0 circr/10]);
p2.DisplayName = '<R>';
p2.LineWidth = 2;
p2.Color =  colour(3,:);
text(txtRng.*cos(txtAngle),txtRng.*sin(txtAngle), names,...
   'HorizontalAlignment','center','VerticalAlignment','bottom',...
   'FontSize', 14, 'FontName', 'Times New Roman')

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
set(findall(gcf, 'String', '  4'),'String', ' Four'); 

ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'normal';
view([-90 -90])