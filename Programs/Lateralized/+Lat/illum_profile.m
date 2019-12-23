function[] = illum_profile(Xfilt, L, R)

%%
[lcolour, rcolour] = colour_palette_leftright(1, 1);

%% check illumination profile

%--- on data ---
%***
figure
<<<<<<< HEAD
polarplot(XLat(:), R(:), 'r.', 'DisplayName', 'right', 'Color', rcolour(2,:))
hold on
polarplot(XLat(:), L(:), 'b.', 'DisplayName', 'left', 'Color', lcolour(2,:))
=======
polarplot(Xfilt(:), R(:), 'r.', 'DisplayName', 'right', 'Color', rcolour(2,:))
hold on
polarplot(Xfilt(:), L(:), 'b.', 'DisplayName', 'left', 'Color', lcolour(2,:))
>>>>>>> master
%polarplot(Xfilt(:), -(L(:)-R(:)), 'k.', 'DisplayName', 'C')
title('on data')
legend

%% illumination profile
% --- theoretical one ----
% --- filled ---
percPm =1;
lum_lin = luminosity_linear(percPm);
theta = (-pi:0.01:pi);
illumL = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),pi-abs(pi-mod(theta,2*pi))); 
illumR = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),abs(pi-mod(theta,2*pi)));

%***
figure
fake = polar(theta, max(illumL)*ones(1,length(illumL)));
set(fake,'Visible','off', 'HandleVisibility', 'off'); 
hold on

<<<<<<< HEAD
pr = polar(theta-pi/2, illumR);
pr.DisplayName= 'right';
pr.LineWidth= 1.5; 
pr.Color = rcolour(3,:);
hold on
xfill = get(pr,'XData');
yfill = get(pr,'YData');
yfillc = yfill;
xfillc = xfill;
yfillc(xfill>0) = [];
xfillc(xfill>0) = [];
yfill(xfill<0) = [];
xfill(xfill<0) = [];
ac = fill(xfillc, yfillc, rcolour(3,:), 'DisplayName', 'conflictual situation (right)');
a = fill(xfill, yfill, rcolour(3,:), 'DisplayName', 'non-conflictual situation (right)');

ac.EdgeColor = rcolour(2,:);
ac.EdgeAlpha = 1;
ac.FaceAlpha =  0.5;

a.EdgeColor = rcolour(2,:);
a.EdgeAlpha = 1;
a.FaceAlpha =  0.2;

hold on
pl = polar(theta-pi/2, illumL);
pl.DisplayName =  'left';
pl.LineWidth = 1.5;
pl.Color =  lcolour(3,:);
hold on
xfill = get(pl,'XData');
yfill = get(pl,'YData');
yfillc = yfill;
xfillc = xfill;
yfillc(xfill>0) = [];
xfillc(xfill>0) = [];
yfill(xfill<0) = [];
xfill(xfill<0) = [];
ac = fill(xfillc, yfillc, lcolour(3,:), 'DisplayName', 'conflictual situation (left)');
a = fill(xfill, yfill, lcolour(3,:), 'DisplayName', 'non-conflictual situation (left)');

ac.EdgeColor = lcolour(2,:);
ac.EdgeAlpha = 1;
ac.FaceAlpha =  0.5;

a.EdgeColor = lcolour(2,:);
a.EdgeAlpha = 1;
a.FaceAlpha =  0.2;
=======
p = polar(theta-pi/2, illumR);
p.DisplayName= 'right';
p.LineWidth= 1.5; 
p.Color = rcolour(3,:);
hold on
a = fill(get(p,'XData'), get(p,'YData'), rcolour(3,:));
a.EdgeColor = rcolour(2,:);
a.EdgeAlpha = 1;
a.FaceAlpha =  0.5;
a.HandleVisibility = 'off';

hold on
p=polar(theta-pi/2, illumL);
p.DisplayName =  'left';
p.LineWidth = 1.5;
p.Color =  lcolour(3,:);
hold on
a = fill(get(p,'XData'), get(p,'YData'), lcolour(3,:));
a.EdgeColor = lcolour(2,:);
a.EdgeAlpha = 1;
a.FaceAlpha =  0.5;
a.HandleVisibility = 'off';
>>>>>>> master

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

<<<<<<< HEAD
ax = gcf;
ax.FontName = 'Times New Roman';
ax.FontSize = 18;
=======
>>>>>>> master

view([-90 -90])

legend
%%

%***
figure
polarplot(theta-pi/2, illumR,...
    'DisplayName', 'right', 'LineWidth', 2, 'Color', rcolour(2,:))

hold on
polarplot(theta-pi/2, illumL,...
    'DisplayName',  'left', 'LineWidth', 2, 'Color',  lcolour(2,:));

legend

set(gca,'ThetaZeroLocation','bottom',...
        'ThetaDir','clockwise')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.ThetaAxisUnits = 'radians';
ax.RAxisLocation = -3*pi/4;    
rticks(0:200:400)
thetaticks((0 : pi/4 : 2*pi))
thetaticklabels({'0', '', '\pi/2', '', '\pi', '', '-\pi/2'})
rtickangle(0)
title(['illumination profile'])

%% lateralized

% --- nice colors for plots ---
blue = [0 0.6 0.7];
red = [0.7 0 0];
orange = [1 0.6 0];
green = [0.2 0.7 0];

% ---
percPm =1;
lum_lin = luminosity_linear(percPm);

theta_deg = 0:360;
theta_rad = deg2rad(theta_deg);

L = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),pi-abs(pi-mod(theta_rad- 3*pi/2,2*pi))); 
R = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),abs(pi-mod(theta_rad- 3*pi/2,2*pi)));

theta_rad_offset = wrapTo2Pi(theta_rad) ;

% ***

fig1 = figure;
%subplot(2,1,1)

plot(theta_rad_offset, L, '-', 'Color', blue, 'Linewidth', 2) 
hold on
plot(theta_rad_offset, R, '-', 'Color', red, 'Linewidth', 2)
plot(theta_rad_offset, L+R, 'Linewidth', 2, 'Color', 0.3*(red+blue))
legend('left intensity', 'right intensity', 'left + right')
ylabel('\muW/cm^2')
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0' '\pi/2' '\pi' '3\pi/2' '2\pi' })

fig2 = figure;
%subplot(2,1,2)

plot(theta_rad_offset, (L-R)/max(L-R), '-k', 'Linewidth', 2)
hold on
area([pi/2 3*pi/2], [1 1], -1,...
    'EdgeAlpha', 0, 'FaceColor', green, 'FaceAlpha', 0.3)
area([0 pi/2], [1 1], -1,...
    'EdgeAlpha', 0, 'FaceColor', orange, 'FaceAlpha', 0.3)
area([3*pi/2 2*pi], [1 1], -1,...
    'EdgeAlpha', 0, 'FaceColor', orange, 'FaceAlpha', 0.3, 'HandleVisibility','off')
plot(pi, 0, 'ksq')
legend('left - right', 'agreement', 'conflict', 'stable point')
ylabel({'relative contrast'; '(I_L-I_R)/I_m_a_x'})
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0' '\pi/2' '\pi' '3\pi/2' '2\pi' })

%% save figure
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/FiguresSerious/';
name1 = 'lateralized_profile1';
name2 = 'lateralized_profile2';
%saveFigurePDF(fig, figpath, name)
saveas(fig1,[figpath name1],'fig')
close(fig1)
saveas(fig2,[figpath name2],'fig')
close(fig2)