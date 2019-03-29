

%%
% --- load spontaneous data ---

timemin = 3;
minBoutNumber = 3;

path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
name = 'PooledSpont.mat';

load([path name], 'Es')


clearvars -except Es
%%
% --- important info ---
pxmm = 11.5;

% --- extract watcha need ---

Xi = Es.Angle;
FishNi = Es.FishN;
Ri = Es.R;
TimeBouti = Es.TimeBout;
xCoordi = Es.Coordinates.x;
yCoordi = Es.Coordinates.y;

% --- modify for watcha need ---

% orientation & trajectories
dXi = diff(Xi,1,2);

dxCoordi = diff(xCoordi, 1, 2);
dyCoordi = diff(yCoordi, 1, 2);
Disti = sqrt(dxCoordi.^2+ dyCoordi.^2);

trajOrientationi = wrapToPi(atan(dyCoordi./dxCoordi));
trajOrientationi(dxCoordi < 0) = trajOrientationi(dxCoordi < 0) - pi;
Alpha = wrapToPi(Xi(:, 1:end-1)) + trajOrientationi;
relAdvancement = cos(Alpha);
Advancementi = Disti.*cos(Alpha);

% time
InterBoutInterval = diff(TimeBouti,1,2);
%%
% --- visualize this shit ---
% chosen sequence
i = 29;

% params for /time and /bout
xl = [10 40];
yl = [5 60];

%***
fig =figure;
visualizeTrajectory_aesthetic...
    (i, Xi, xCoordi/pxmm, yCoordi/pxmm, Disti/pxmm, Advancementi/pxmm, dXi, trajOrientationi)

xlim(xl)
ylim(yl)
xlabel('mm')
ylabel('mm')

%% save figure
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Figures201808/';
name = ['spontaneous_trajectory_xytheta'];
%saveFigurePDF(fig, figpath, name)
saveas(fig,[figpath name],'fig')
%% 
% --- in time ---
loadplease = 1;
if loadplease
    if 0
        [Dates, Fish, exptype, ~, intmax] = chooseExpType('exp60');
        p = 2; q = 1;
        fish = char(Fish(p,q));
        date = char(Dates(p));
        
        load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype '/' intmax '/' date '/' fish '/data' fish '.mat']);
        
        [angle_ini, coordinates_ini, finfo] = remove_sequence_ini(D, 3);
        [angleCum_ini] = angle_cum(angle_ini, finfo);
    end
    i = 11;
    framerate_ini = finfo(i,3);
    theta = angle_ini(i,:);
    theta_cum = angleCum_ini(i,:);
    theta_cum(isnan(theta_cum)) = [];
    theta_cum_sm = smooth(theta_cum, 20, 'sgolay');
    dth = diff(theta_cum_sm);
    
    xtemp = coordinates_ini(:,1,i);
    ytemp = coordinates_ini(:,2,i);
    xtemp=squeeze(xtemp);
    ytemp=squeeze(ytemp);
    xtemp(xtemp==0) = nan;
    ytemp(ytemp==0) = nan;
    xtemp(isnan(xtemp)) = [];
    ytemp(isnan(ytemp)) = [];
    xtemp_sm = smooth(xtemp, 20, 'sgolay');
    ytemp_sm = smooth(ytemp, 20, 'sgolay');
    dx = diff(xtemp_sm);
    dy = diff(ytemp_sm);
end
dd = sqrt(dx.^2+dy.^2);

unitvecty1 = 1*sin(deg2rad(theta_cum_sm)-pi);
unitvectx1 = 1*cos(deg2rad(theta_cum_sm));
%% *** without velocities
fig = figure;
v = VideoWriter([figpath 'spontaneous_trajectory.avi']);
v.FrameRate = 74;
v.Quality = 100;
open(v);
for i = 1 : length(xtemp)
    plot(xtemp(1:i)/pxmm, ytemp(1:i)/pxmm, 'Color', [0 0.2 0.2])
    hold on
    uvx1 = linspace(xtemp(i)/pxmm, xtemp(i)/pxmm+unitvectx1(i), 10);
    uvy1 = linspace(ytemp(i)/pxmm, ytemp(i)/pxmm+unitvecty1(i), 10);
    drawArrow([uvx1(1) uvx1(end)], [uvy1(1) uvy1(end)]);
    text(35, 57, [num2str(floor(10*i/framerate_ini)/10) ' secs'])
    set(gca,'DataAspectRatio',[1,1,1])
    xlabel('mm')
    ylabel('mm')
    xlim(xl)
    ylim(yl)
    frame = getframe(fig);
    writeVideo(v,frame);
    hold off
    %pause(1/framerate_ini/10)
end
close(v)

%% *** with velocities
figure
set(gca,'DataAspectRatio',[1,1,1])
xlabel('mm')
ylabel('mm')
% v = VideoWriter([figpath 'StereoDistributionsOverTime.avi']);
%      v.FrameRate = 2;
% open(v);
for i = 1 : length(xtemp)
    subplot(2,1,1)
    plot(xtemp(1:i)/pxmm, ytemp(1:i)/pxmm, 'k')
    text(35, 57, [num2str(floor(10*i/framerate_ini)/10) ' secs'])
    xlim(xl)
    ylim(yl)
    %frame = getframe(fig);
    %writeVideo(v,frame);
    
    subplot(2,1,2)
    yyaxis left
    plot((1:i)/framerate_ini, dd(1:i)/pxmm*framerate_ini)
    ylim([0 20])
    
    yyaxis right
    plot((1:i)/framerate_ini, deg2rad(dth(1:i))*framerate_ini);
    ylim([-pi pi])
    
    if i/framerate_ini < 100/framerate_ini
        xlim([0 200/framerate_ini])
    else
        xlim([(i-100)/framerate_ini (i+100)/framerate_ini])
    end
    
    pause(1/framerate_ini)
end
