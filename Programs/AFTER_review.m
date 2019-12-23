
%% select trajectories in simulated data
load('/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/AnalysisOutput/trajectories/no_stim.mat')
sel = [3 15 64 67 68 94 99 100 109 117 124 127 143];

%% select trajectories in spontaneous data
load('/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Annexe/ForGuillaume/SpontaneousSwimCoordinates.mat')
xdata = Coordinates.x;
ydata = Coordinates.y;
seldata = [2 6 29 26 32 40 45 51 55 57]; %12

xdatasel = xdata(seldata,:);
ydatasel = ydata(seldata,:);

xdatasel = xdatasel - xdatasel(:,1);
ydatasel = ydatasel - ydatasel(:,1);

seq_lengths = nan(size(xdatasel,1),1);
for i = 1 : size(xdatasel,1)
    seq_lengths(i) = find(isnan(xdatasel(i,:)),1);
end
xdatasel = xdatasel(:, 1:max(seq_lengths));
ydatasel = ydatasel(:, 1:max(seq_lengths));

subplot(2,1,1)
plot((xdatasel/11.5)', (ydatasel/11.5)',...
    '.-', 'Linewidth', 1.2, 'Markersize', 8, 'MarkerFaceColor', [0 0 0])
hold on
plot(0, 0, 'pentagram', 'Markersize', 12,...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0])
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Times New Roman';
xticks([-50 : 50 : 50])
yticks([-50 : 50 : 50])
ylabel('[mm]')
xlabel('[mm]')
axis equal
title('Data')

%% seb's trajectories (load figure of trajectory with IBI = 1sec)
ax=gca;
%xseb = ax.Children(2).XData;
%yseb = ax.Children(2).YData;
hold off
subplot(2,1,2)
for i = 1 : length(seq_lengths)
    start_point = round(rand*length(xseb));
    end_point = start_point + seq_lengths(i);
    seq_lengths(i)
    xdisp = xseb(start_point:end_point)*10;
    xdisp = xdisp-xdisp(1);
    ydisp = yseb(start_point:end_point)*10;
    ydisp = ydisp-ydisp(1);
    plot(xdisp, ydisp, ...
        '.-', 'Linewidth', 1.2, 'Markersize', 8, 'MarkerFaceColor', [0 0 0])
    hold on
end
hold on
plot(0, 0, 'pentagram', 'Markersize', 12,...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0])
title('Simulation')
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Times New Roman';
xticks([-50 : 50 : 50])
yticks([-50 : 50 : 50])
ylabel('[mm]')
xlabel('[mm]')
axis equal
