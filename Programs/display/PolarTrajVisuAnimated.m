%% format and save data for Raphael's visualization

%% lateralized data
%% --- load whole dataset ---
Lat.a.Load
time = TimeBout;
alpha = XLat-pi/2;
x = xCoord;
y = yCoord;
l = Dist/pxmm;

%% --- or individual experiment ---
% load the experiment
pxmm = 11.5;

load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
    'lateralisation/18-10-05/fish1.mat'], 'E')
alpha1 = deg2rad(E.experiment.angleCum);
x1 = squeeze(E.experiment.coordinates(:,1,:))'/pxmm;
y1 = squeeze(E.experiment.coordinates(:,2,:))'/pxmm;
%l1 = sqrt(diff(x1,1,2).^2+diff(y1,1,2).^2);
l1 = sqrt( (x1-repmat(x1(:,1),[1 size(x1,2)])).^2 + (y1-repmat(y1(:,1),[1 size(y1,2)])).^2 );
fr1 = E.experiment.framerate(:,3);
ts1 = ((1./fr1).*[0:size(alpha1,2)-1]);

load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
    'lateralisation/18-09-28/fish3.mat'], 'E')
alpha2 = deg2rad(E.experiment.angleCum);
x2 = squeeze(E.experiment.coordinates(:,1,:))'/pxmm;
y2 = squeeze(E.experiment.coordinates(:,2,:))'/pxmm;
%l2 = sqrt(diff(x2,1,2).^2+diff(y2,1,2).^2);
l2 = sqrt( (x2-repmat(x2(:,1),[1 size(x2,2)])).^2 + (y2-repmat(y2(:,1),[1 size(y2,2)])).^2 );
fr2 = E.experiment.framerate(:,3);
ts2 = ((1./fr2).*[0:size(alpha2,2)-1]);

how_many_nans_to_add = abs(size(alpha2,2)-size(alpha1,2));
nans_to_add = NaN(size(alpha1,1), how_many_nans_to_add);

alpha = [alpha1 nans_to_add; alpha2];
l = [l1 nans_to_add; l2];
ts = [ts1 nans_to_add; ts2];
% rho = cumsum(l,2);
% rho = [zeros(size(rho,1),1) rho];
rho = l;

arenans = isnan(alpha)+1;
arenans(arenans==2) = NaN;
ts = ts.*arenans;
alpha = alpha-3*pi/2;

%% correct alpha
[too_bigr,too_bigc]  = find(abs(diff(alpha,1,2)) > pi/2);
for i = 1 : length(too_bigr)
    alpha(too_bigr(i), too_bigc(i)+1) = (alpha(too_bigr(i), too_bigc(i)-1) + alpha(too_bigr(i), too_bigc(i)+3))/2;
end
%% --- visualize ---

alphasmooth = alpha;
for i = 1 : size(alpha,1)
    alphasmooth(i,:) = smooth(alpha(i,:), 7, 'moving');
end


% --- metric = time ---
for i = 1 : 1000 %size(rho,2)
    for j = 1 : size(rho,1)
        polarplot(alpha(j,1:i), rho(j,1:i))
        hold on
    end
    drawnow
%     if mod(i,10)==0
%         disp(nanmean(ts(:,i)))
%     end
    hold off
end

% --- metric = bout ---
% for i = 1 : size(alpha,2)
%     polarplot(alpha(:,1:i), rho(:,1:i))
%     drawnow
% end

%% save
trajectories.alpha.comment = 'orientation towards source (rad)';
trajectories.alpha.values = alpha;
trajectories.alpha.smoothedvalues = alphasmooth;
trajectories.rho.comment = 'distance from origin (mm)';
trajectories.rho.values = rho;
trajectories.time.comment = 'time stamp (sec)';
trajectories.time.values = ts;

save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/TrajectoryDisplay/trajectories.mat'],...
    'trajectories')