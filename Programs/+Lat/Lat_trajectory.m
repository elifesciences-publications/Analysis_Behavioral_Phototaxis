
% use after AnalysisLat

i = 24;
fish = find(FishN == i);
x = XLat(fish,:);
t = TimeBout(fish,:);
%***
figure;
polarhistogram(x(:), 30)

%***
figure
visualizeTrajectoryAndAdv(fish, Xlab, xCoord, yCoord, Dist, Advancement, dX, trajOrientation, 0)

meanX = circ_mean(x);
R = circ_r(x);
Rproj = R.*cos(meanX+pi);

%***
figure;
plot(Rproj(1:30))

%%
size(x,1)
numel(x) - sum(isnan(x(:)))

%% ***
% r = linspace(0, 0.8, size(x,1));
% g = linspace(0, 0.8, size(x,1));
% b = linspace(0, 0.8, size(x,1));
 colour = [r' g' b'];
 orange = [1 0.6 0];

fig = figure;
%p = plot(t', x', 'k'); % time
plot([0 35], [pi pi], '--', 'Color', orange, 'Linewidth', 1.5)
hold on
plot([0 35], [-pi -pi], '--', 'Color', orange, 'Linewidth', 1.5)
plot(x'-pi/2, 'k'); % bout #
xlim([1 15])
xlabel('bout #')
ylabel('cumulative orientation before bout (rad)')
ylim([-3*pi 3*pi])
yticks([-2*pi : pi : 2*pi])
yticklabels({'-2\pi' '-\pi' '0' '\pi' '2\pi'})

%% save figure
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/FiguresSerious/';
name = ['lateralized_trajectoriesVSbouts'];
%saveFigurePDF(fig, figpath, name)
saveas(fig,[figpath name],'fig')

%% temporal
Lat_DatesFish
p=18; q=1;
fish = char(Fish(p,q));
date = char(Dates(p));
load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
    exptype '/' date '/' fish '/data' fish '.mat']);

acum = D.experiment.angleFiltered;
acum(acum == 0) = NaN;

plot(acum')

