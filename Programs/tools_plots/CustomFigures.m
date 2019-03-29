% custom plots for presentation

%% tropotaxis experiment
x = 0 : 360;
righteye = [linspace(1, 0, floor(length(x)/2))...
    linspace(0, 1, ceil(length(x)/2))];
lefteye = [linspace(0, 1, floor(length(x)/2))...
    linspace(1, 0, ceil(length(x)/2))];
fig = figure;
ax1 = plot(deg2rad(x), righteye);
hold on
ax2 = plot(deg2rad(x), lefteye);
ax3 = plot(deg2rad(x), righteye+lefteye);

yticks([0 0.5 1])
ylim([0 1.05])

xlim([0 deg2rad(x(end))])
xticks([0:pi/2:2*pi])
xticklabels({'0' '\pi/2' '\pi' '3\pi/2' '2\pi'})

xlabel('orientation (rad)')
ylabel('normalized intensity')
l = legend('left', 'right', 'left+right');

fig. Color = [1 1 1];
ax1.LineWidth = 2;
ax2.LineWidth = 2;
ax3.LineWidth = 2;
l.FontSize = 16;

%% klinotaxis experiment
lumRad = deg2rad(lum_th(:,1));
lumW = lum_th(:,2);

%***
f = figure;
f.Name = 'Intensity profile';

subplot1 = subplot(1,2,1);
plot(lumRad, 1000*lumW, 'Linewidth', 2);
xlabel('orientation (rad)', 'FontSize', 16);
ylabel('intensity (mW/cm^2)', 'FontSize', 16);
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlim([0 pi])
title('')

subplot2 = subplot(1,2,2);
plot(1000*lumW(2:end), diff(1000*lumW), 'Linewidth', 2);
xlabel('intensity (mW/cm^2)', 'FontSize', 16);
ylabel('\delta instensity (mW/cm^2)', 'FontSize', 16);
xlim([0 1000*lumW(end)])
title('\Delta I/I constant')

box(subplot1,'off');
box(subplot2,'off');
set(subplot2,'FontSize',14,'YTick',[0 : 4]);
set(subplot1,'FontSize',14);

% subplot(2,2,3:4)
% plot(lumRad(2:end), diff(1000*lumW)./lumW(2:end)*1000)
% xlabel('orientation (rad)');
% xticks([0 pi/4 pi/2 3*pi/4 pi])
% xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
% xlim([0 pi])
% ylabel('dI/I');

%% trajectory in time
lumiere = lum_th(:,3);
luorient = deg2rad(lum_th(:,1));

ex_traj = deg2rad(angleSourceCorr(4,:));
ex_coordx = coordinates(:,1,4);
ex_coordx(ex_coordx==0)=NaN;
ex_coordy = coordinates(:,2,4);
ex_coordy(ex_coordy==0)=NaN;

fig = figure;

% --- plot gradient ---
subplot(1,2,1);
gradrect = pi/64;
for i = 0 : gradrect : max(ex_traj) - gradrect
    c = -(cos(i)+1)/8+1;%lumiere(rad2deg(i));
    colorgrad = [c c c];
    rect_H = rectangle('Position', [0, i, l, gradrect]); 
    set(rect_H, 'FaceColor', colorgrad);
    set(rect_H, 'EdgeColor', colorgrad) 
end

hold on

w = 100;
l = length(ex_traj)-sum(isnan(ex_traj));

for i = 1 : w : l -w
    c = [0.3 0.3 (i/l)/1.5];
    
    subplot(1,2,1);
    plot(i:i+w, ex_traj(i:i+w), 'Color', c);
    hold on
    
    subplot(1,2,2);
    plot(ex_coordx(i:i+w), ex_coordy(i:i+w), 'Color', c, 'Linewidth', 2)
    hold on
end

subplot1 = subplot(1,2,1);
set(subplot1,'FontSize',14)
xlabel('frames')
xlim([0 l])
ylim([0 max(ex_traj)])
ylabel('cumulative orientation relative to the source')
yticks([0 : pi : 10*pi])
yticklabels({'0' '\pi' '2\pi' '3\pi'})

subplot2 = subplot(1,2,2);
xlabel('x-coordinates (px)')
ylabel('y-coordinates (px)')
set(subplot2,'FontSize',14)
%%
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,'trajectory_biased_exp60','-dpdf','-r0')

%% find bouts
BoutSpot(coordinates, angleSourceCorr, framerateInfo, 1)
%% trajectory per bout
ex_traj_bouts = deg2rad(angleBoutSource(4,:));

fig =1;
BoutSpot(coordinates, angleLabCumCorr_ini, framerateInfo_ini, fig);

xlabel('time (seconds)')
xlim([0 323])
ylabel('orientation relative to the source')
yticks([0:2*pi:10*pi])
yticklabels({'0' '2\pi' '4\pi' '6\pi' '8\pi' '10\pi'})

legend('recorded orientation', 'smoothed orientation', 'bouts')

fig = figure;
plot(ex_traj_bouts, 'Linewidth', 2, 'Color', colorgris)
ylabel('cumulative orientation relative to the source')
yticks([0 : 2*pi : 10*pi])
yticklabels({'0' '2\pi' '4\pi' '6\pi' '8\pi' '10\pi'})
xlabel('bout #')
ax =gca;
ax.FontSize = 14;

%%
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,'bouttrajectory_biased_exp60','-dpdf','-r0')
