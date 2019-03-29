clc

% ===¨Parameters ==========================================================
% from data
%TrajFile = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/TrajectoryDisplay/trajectories.mat';
% from simu
TrajFile = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Programs/TrajectoryDisplay/trajectories.mat';

OutDir = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/TrajectoryDisplay/Frames';

speed = 1;    % Speed factor
dur = 60;        % Movie duration, in experimental time (seconds)
pers = 20;     % Persistence time, in experimental time (seconds)

% =========================================================================

% --- Data ----------------------------------------------------------------

% --- Load
FMat = load(TrajFile);
T = FMat.trajectories.time.values;
Rho = FMat.trajectories.rho.values;
Alpha = FMat.trajectories.alpha.values;
n = size(T,1);

% --- Regularization of vector sizes
tmp = min([size(T,2) size(Rho,2) size(Alpha,2)]);
T = T(:,1:tmp);
Rho = Rho(:,1:tmp);
Alpha = Alpha(:,1:tmp);
rhomax = max(Rho(:));

% --- Interpolation

% Number fo frames
Nf = dur*25/speed + 1;

t = linspace(0, dur, Nf);
% rho = NaN(n, Nf);
% alpha = NaN(n, Nf);
% for i = 1:n
%     K = ~isnan(Rho(i,:));
%     rho(i,:) = interp1(T(i,K), Rho(i,K), t);
%     alpha(i,:) = interp1(T(i,K), Alpha(i,K), t);
% end
rho = Rho;
alpha = Alpha;

% --- Display loop --------------------------------------------------------

% --- Preparation
cm = rand(n,3);
if ~exist(OutDir, 'dir')
    mkdir(OutDir);
end

set(gcf, 'color', 'white');

for k = 1:Nf
    
    t = (k-1)/25/speed;
 
    clf
    polaraxes('ThetaZeroLocation', 'top', 'RLim', [0 rhomax]);
    hold on
    
    I = max(1,round(k-pers*25/speed)):k;
    for i = 1:n
        polarplot(alpha(i,I), rho(i,I), 'color', cm(i,:))
        polarscatter(alpha(i,k), rho(i,k), ...
            'MarkerFaceColor', cm(i,:), ...
            'MarkerEdgeColor', 'none')
    end
    
    title(['t = ' num2str(t, '%.02f') 's']);
    
    F = getframe(gcf);
    imwrite(F.cdata, [OutDir filesep 'Frame_' num2str(k, '%06i') '.png']);
    
end
