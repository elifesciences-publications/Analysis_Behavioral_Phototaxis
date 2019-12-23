%% TEMPORAL PHOTOTAXIS DATA
% different experiment types separetely

%% select experiments

%experiment = 'exp60';
experiment = 'exp30';
%experiment = 'sin60';

[Dates, Fish, exptype, lum_th, intmax] = Temp.perexperiment.select_exp_type(experiment);

%% ************************************************************************
...FETCH AND POOL ALL INDIVIDUAL EXPERIMENTS...

% --- some parameters to tune ---
timemin = 3; 
minBoutNumber = 3;

    % initial LOOP
[L,l] = initial_loop(Dates, Fish, exptype, intmax, timemin);

    % pooling LOOP
[Xrad, Xlab, Lu, TimeBout, xCoord, yCoord, FishN] = ...
    Temp.pooling.loop(Dates, Fish, l, L, exptype, intmax, timemin, minBoutNumber);


if ( [L, l] - size(Xrad)) < 0
    warning('problem with array size initial estimation')
end
%% ************************************************************************
... Clean up data, create useful variables

% --- wrap & create dX
Xwrapped_0pi = abs( wrapToPi ( Xrad ) ); % wrap to [0 pi]

Xwrappeds0 = pi - Xrad ;
Xwrappeds0 = abs( wrapToPi( Xwrappeds0 ) );

dX = diff( Xrad, 1, 2 );
dXw = diff( Xwrapped_0pi, 1, 2 );

XforDiff = Xrad(:,1:end-1);
XwrappedForDiff = abs(wrapToPi(Xrad(:,1:end-1)));

dLu = diff(Lu,1,2);
LuForDiff = Lu(:,1:end-1);

% luminosity
lumRad = deg2rad(lum_th(:,1));
lumW = lum_th(:,2);

% coordinates
dxCoord = diff(xCoord, 1, 2);
dyCoord = diff(yCoord, 1, 2);
Dist = sqrt(dxCoord.^2+ dyCoord.^2);

%advancement
trajOrientation = wrapToPi(atan(dyCoord./dxCoord));
trajOrientation(dxCoord < 0) = trajOrientation(dxCoord < 0) - pi;
dAlpha = wrapToPi(Xlab(:, 1:end-1)) + trajOrientation;
R = Dist.*cos(dAlpha);

% some stats on sequences
    % number of fish
seqPerFish = histc(FishN, unique(FishN));
NFish = length(seqPerFish);
meanNofSeqPerFish = mean(seqPerFish);

NofSeqs = size( Xrad, 1 );
boutsPerSeq = size( Xrad, 2 ) - sum(isnan(Xrad),2);
maxnseq = max(seqPerFish);
%% ************************************************************************

%% VISUALIZE TRAJECTORIES with ADVANCEMENT R

% --- look at individual trajectories ---
fish = 10:20;

visualizeTrajectoryAndAdv(fish, Xlab, xCoord, yCoord, Dist, R, dX, trajOrientation)

%% illumination profile

%***
fig = figure;
fig.Name = 'Intensity profile';

s = subplot(2,2,1);
plot(lumRad, lumW);
xlabel('fish orientation (rad)');
ylabel('intensity (W)');
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlim([0 pi])
title('')

subplot(2,2,2)
plot(lumW(2:end), diff(lumW));
xlabel('intensity (W)');
ylabel('\delta instensity (W/rad)');
xlim([0 lumW(end)])

subplot(2,2,3:4)
plot(lumRad(2:end), diff(lumW)./lumW(2:end))
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlabel('orientation (rad)');
ylabel('\deltaI/I (rad^-^1)');
xlim([0 pi])

title(s, experiment)

%%
name = ['all_illum'];
path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Figures/';
saveFigurePDF(fig, path, name)

%% distribution of X
colorgris = [22/256 73/256 132/256];

%***
fig = figure;
subplot(2,1,1);
polarhistogram(wrapToPi(Xrad(:)),0 : pi/12 :  2*pi,...
    'Normalization', 'probability', 'FaceColor', colorgris,'FaceAlpha',.3);
set(gca,'ThetaZeroLocation','left',...
        'ThetaDir','clockwise');
thetaticks( 0 :90: 360 )
thetaticklabels({'0' '\pi/2' '\pi' '3\pi/4'})
rticks([0.03 0.06])
title('distribution of X at each bout, X \in [-\pi \pi]')
circr = circ_r(wrapToPi(Xrad(:)));
circmean = circ_mean(wrapToPi(Xrad(:)));
hold on
polarplot([0 circmean], [0 circr], 'Linewidth', 2, 'Color', [0.8 0.3 0.3]) 
ax = gca;
ax.RAxisLocation = 45;
ax.FontSize = 14;
ax.FontName = 'Helvetica Neue';
ax.FontWeight = 'normal';

subplot(2,1,2)
histogram( Xwrapped_0pi, 0 : pi/32 :  pi , 'FaceColor', colorgris,...
    'Normalization', 'probability')
ax = gca;
ax.FontSize = 14;
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
yticks(0:0.01:0.05)
xlim([0 pi])
title('distribution of X, X\in [0 \pi]')
ax.FontName = 'Helvetica Neue';
%%
name = [experiment '_Xdistribution_biased'];
saveFigurePDF(fig, path, name)

%% misc stats ORIENTATION

% --- mean and variance of dX for different orientations ---
w = pi/16 ;
bins = -pi : w : pi ;
meandX = NaN(length(bins), 1);
vardX = NaN(length(bins), 1);
meandist = NaN(length(bins), 1);
semdX = NaN(length(bins), 1);
stddX =  NaN(length(bins), 1);
n = NaN(length(bins), 1);
nL =  NaN(length(bins), 1);
nR =  NaN(length(bins), 1);
for i = 2 : length(bins)
    dXselected = dX(wrapToPi(XforDiff) < bins(i) & wrapToPi(XforDiff) >= bins(i-1));
    meandX(i-1) = nanmean((dXselected));
    vardX(i-1) = nanvar(dXselected);
    stddX(i-1) = nanstd(dXselected);
    n(i-1) = length(dXselected);
    nL(i-1) = sum(dXselected>0);
    nR(i-1) = sum(dXselected<0);
    semdX(i-1) = nanstd(dXselected)/sqrt(n(i-1));
    %figure;
    %histogram(dXselected, [-3:0.1:3])
%      figure
%      histogram(dXselected)
%      drawnow
%      pause(1)
end

%***
fig = figure;

%scatter(XwrappedForDiff(:), abs(dX(:)), 12)
subplot(2,1,1)
grid on
plot(bins + w/2, meandX, bins + w/2, vardX, 'Linewidth', 2)
hold on
errorbar(bins + w/2, meandX, semdX, '.', 'Color', [0.2 0.2 0.2])
xlabel('orientation (rad)')
ylabel('angular displacement per bout dX')
legend( '<(\delta X)>', '<\delta X ^2>')
xticks([-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'-\pi', '-3\pi/4', '-\pi/2', '-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
xlim([-pi pi])

ax = gca;
ax.FontSize = 14;

subplot(2,1,2)
grid on
plot(bins + w/2, n, 'Linewidth', 2)
hold on
plot(bins + w/2, nL, 'r', bins + w/2, nR, 'b')
plot(bins + w/2, nL-nR, 'k')
ylabel('bout number')
legend('all bouts', 'left bouts', 'right bouts', 'L-R')

ax = gca;
ax.FontSize = 14;

%% save figure
name = [experiment '_A_D_vs_orient'];
saveFigurePDF(fig, path, name)

%% ---  distributions of dX in broad bins ---
w = pi/4 ; % bin size
bins = -pi : w : pi ;
bins = 0 : w : pi;

%***
figure
for i = 2 : length(bins)
    dXselected = dX(abs(wrapToPi(XforDiff)) < bins(i) & abs(wrapToPi(XforDiff)) >= bins(i-1));
    if isempty(dXselected)
        xdens = []; dens = [];
    else
        [dens, xdens] = ksdensity( dXselected,...
        'Function', 'pdf', 'Kernel', 'Normal'); % , 'Support', 'positive'
    end
    plot(xdens, dens, 'Linewidth', 2)
    hold on
end
%xlim([0 max(dXi(:))])
xlabel('dX')
ylabel('distribution')
legend(num2str(bins'))
%% --- luminosity ---
bins = linspace(lumW(1), lumW(end), 24);
bins = linspace(min(dLu(:)), max(dLu(:)), 30);
w = mean(diff(bins));

meandXl = NaN(length(bins), length(unique(FishN)));
vardXl = NaN(length(bins), length(unique(FishN)));
semdXl = NaN(length(bins), length(unique(FishN)));

for j = 1 : length(unique(FishN))
    
    for i = 2 : length(bins)
        dXfish = dX(FishN == j, :);
        %lufish = LuForDiff(FishN == j, :);
        lufish = dLu(FishN == j, :);
        dXselected = dXfish(lufish < bins(i) & lufish> bins(i-1));
        meandXl(i-1,j) = nanmean(dXselected);
        vardXl(i-1,j) = nanvar(dXselected);
        n = length(dXselected);
        semdXl(i-1,j) = nanstd(dXselected)/sqrt(n);
        %figure
        %histogram(dXselected)
    end
    
end
%***
fig = figure;
%scatter(LuForDiff(:), abs(dX(:)), 12)

subplot(1,2,1)
plot(bins + w/2, nanmean(meandXl'), 'Linewidth', 2)
hold on
errorbar(bins + w/2, nanmean(meandXl'), nanstd(meandXl')/length(unique(FishN)), '.', 'Color', [0.2 0.2 0.2])
legend('mean dX', 'sem on n. of fish')
ax = fig.CurrentAxes;
ax.YGrid = 'on';
ax = gca;
ax.FontSize = 14;
xlabel('lum (W)')
ylabel('dX')
title('mean')

subplot(1,2,2)
plot(bins + w/2, nanmean(vardXl'), 'Linewidth', 2)
hold on
errorbar(bins + w/2, nanmean(vardXl'), nanstd(vardXl')/length(unique(FishN)), '.', 'Color', [0.2 0.2 0.2])
legend('var dX', 'sem on n. of fish')
title('var')
xlabel('lum (W)')
ylabel('dX^2')
ax = fig.CurrentAxes;
ax.YGrid = 'on';
ax = gca;
ax.FontSize = 14;

%% --- distributions of dX in broad luminosity bins ---
bins = linspace(lumW(1), lumW(end), 10);
w = mean(diff(bins));

%***
figure
for i = 2 : length(bins)
    dXselected = dX(LuForDiff < bins(i) & LuForDiff > bins(i-1));
    if isempty(dXselected)
        xdens = []; dens = [];
    else
        [dens, xdens] = ksdensity( dXselected,...
        'Function', 'pdf', 'Kernel', 'normal', 'bandwidth', 0.1); % , 'Support', 'positive'
        size(dens)

    end
    if i==2 || i == 4 || i == 8
        plot(xdens, dens, 'Linewidth', 2)
        hold on
    end
end
%xlim([0 max(dXi(:))])
xlabel('dX')
ylabel('pdf')
title('distrib pour differentes luminosit?s')
legend(num2str(bins([1, 2])),...
    num2str(bins([3, 4])),...
    num2str(bins([7, 8])) )


%% Bias : A
Var = dXw;
opa = 1;

Aall = nanmean(Var, 2);
% --- A per fish & per sequence ---

Afish = NaN(length(seqPerFish), maxnseq);
boutsPerFish = NaN(length(seqPerFish), 1);
Rfish = NaN(length(seqPerFish), 1);
CircMeanfish = NaN(length(seqPerFish), 1);
s=1;
check=[];
for fish = 1 : length(seqPerFish)
    rowsOfsequencesOneFish = s : s + seqPerFish(fish) - 1;
    afish = nanmean(Var(rowsOfsequencesOneFish, :), 2);
    Afish(fish, 1 : length(rowsOfsequencesOneFish) ) = afish;
    dXOnefish = Var(rowsOfsequencesOneFish,:);
    boutsPerFish(fish) = numel(dXOnefish) - sum(isnan(dXOnefish(:)));
    linXfish = LinearizeSeqForOneIndiv(Xrad, rowsOfsequencesOneFish);
    Rfish(fish) = circ_r(linXfish);
    CircMeanfish(fish) = circ_mean(linXfish);
    s = s + seqPerFish(fish); 
end

% --- mean A per fish with sem ---
meanAfish = nanmean(Afish, 2);
semAfish = nanstd(Afish,1,2)./sqrt(boutsPerFish);

%fig = figure;
%fig.Name = 'Drift during biased swimming';

subplot(1,2,1)
bar(meanAfish, 'FaceColor', colorgris, 'FaceAlpha', opa);
hold on
errorbar(meanAfish, semAfish, '.', 'Color', [0.7 0 0.2])
xlabel('fish #')
ylabel('A = <\delta X > (rad)')
xlim([1 14])
ylim([-0.12 0.12])
legend('mean per fish', 'sem')
ax = gca;
ax.FontSize = 14;

subplot(1,2,2)
histogram(Aall, [-1:0.05:1],...
    'Normalization', 'probability', 'FaceColor', colorgris, 'FaceAlpha', opa)
text(5, 30, ['Mean = ' num2str(mean(Aall)) ' rad'])
xlabel('radians')
if opa == 1
    text(0.2, 0.25, [{'Mean dX wrapped = '} {num2str(nanmean(Aall)) ' rad'}])
else
    text(0.2, 0.12, [{'Mean dX = '} {num2str(nanmean(Aall)) ' rad'}])
end
ax = gca;
ax.FontSize = 14;


%% Autocorrelations of dX
[ACnorm] = xcorrMatrixRows (dX);

% ***
figure
plot(0:length(ACnorm)-1, ACnorm)
xlim([0 50])
grid on
title('Biased - autocorrelation of dX')

%% Diffusion : D

% --- D per fish ---
% s=1;
% for i = 1 : length(seqPerFish)
%     figure
%     [MSDps] = MSDperseq( Xrad(s : s + seqPerFish(i) - 1 , :));
%     plot(MSDps')
%     s = s + seqPerFish(i); 
% end

% --- D for all sequences ---

bias = Aall;

X_X0 = Xrad - Xrad(:,1);
sd0 = (X_X0 - bias).^2;
msd0 = mean(sd0, 1, 'omitnan');

[MSD] = MSDXt0shuffled ( Xrad , 0);

% ***
figure
plot(MSD);
ll = 5;
ul = 60;
title('Biased - < (X(t)-X(t0)) ^2>')
hold on
p = polyfit( (ll:ul), MSD(ll:ul) ,1); 
%plot(1:ul, (1:ul)*p(1)+p(2))
text(5, 20, ['pente = ' num2str(p(1))])

%% advancement

pxmm = 11.5;

% ***
fig = figure;
subplot(3,1,1)
plot(dX(:), R(:)/pxmm, '.', 'MarkerSize', 0.1)
xlabel('dX')
ylabel('advancement (mm)')
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi', '-\pi/2','0','\pi/2','\pi'})
xlim([-pi pi])
subplot(3,1,2)
plot(wrapToPi(XforDiff(:)), R(:)/pxmm, '.', 'MarkerSize', 0.1)
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi', '-\pi/2','0','\pi/2','\pi'})
xlim([-pi pi])
xlabel('orientation')
ylabel('advancement (mm)')
subplot(3,1,3)
plot(LuForDiff(:), R(:)/pxmm, '.', 'MarkerSize', 0.1)
xlim([0 0.22])
xlabel('luminosity')
ylabel('advancement (mm)')

%%
name = [experiment '_Rscatter'];
saveFigurePDF(fig, path, name)
%%
w = pi/16 ;
bins = 0 : w : pi ;
meanR = NaN(length(bins), 1);
varR = NaN(length(bins), 1);
semR = NaN(length(bins), 1);
for i = 2 : length(bins)
    Rselected = R(abs(wrapToPi(XforDiff)) < bins(i) & wrapToPi(XforDiff) >= bins(i-1))/pxmm;
    meanR(i-1) = nanmean(Rselected);
    varR(i-1) = nanvar(Rselected);
    n = length(Rselected);
    semR(i-1) = nanstd(Rselected)/sqrt(n);
    %figure;
    %histogram(dXselected, [-3:0.1:3])
end

%***
fig = figure;
%scatter(abs(wrapToPi(XforDiff(:))), R(:), 12)
hold on
plot(bins + w/2, meanR, 'Linewidth', 2)
errorbar(bins + w/2, meanR, semR, '.', 'Color', [0.2 0.2 0.2])
xlabel('orientation (rad)')
ylabel('advancement per bout (mm)')
%legend( 'indiv', 'mean R', 'sem')
legend('mean R', 'sem')
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlim([0 pi])
ax = gca;
ax.FontSize = 14;
grid on
%%
name = [experiment '_RmeanOrient'];
saveFigurePDF(fig, path, name)
%%
% advancement VS luminosity
bins = linspace(0, max(Lu(:)), 50);
w = mean(diff(bins));

meandRla = NaN(length(bins), 1);
vardRla = NaN(length(bins), 1);
semdRla = NaN(length(bins), 1);
for i = 2 : length(bins)
    Rselected = R(LuForDiff < bins(i) & LuForDiff> bins(i-1))/pxmm;
    n = length(Rselected);
    semdRla(i-1) = nanstd(Rselected)/sqrt(n);
    meandRla(i-1) = nanmean(Rselected);
    vardRla(i-1) = nanvar(Rselected);
end

%***
figure
plot(bins + w/2, meandRla, 'Linewidth', 2);%, bins + w/2, vardRla, 'Linewidth', 2)
hold on
errorbar(bins + w/2, meandRla, semdRla, '.', 'Color', [0.2 0.2 0.2])
legend( '<\delta X>', 'sem')
ylabel('advancement per bout (mm)')
xlabel('lum')
ax = gca;
ax.YGrid = 'on';
ax.FontSize = 14;

%%
name = [experiment '_RmeanLum'];
saveFigurePDF(fig, path, name)
%% memory : --- single bout ---

Vart1 = dLu(:, 1:end-1);
Vart2 = dX(:, 2:end).*(dX(:, 1:end-1)./abs(dX(:, 1:end-1)));
%Vart1 = shuffle(Vart1);


b = max( abs(min(Vart1(:))), max(Vart1(:)));
bins = linspace(-b, b, 50);
w = mean(diff(bins));

meanV2 = NaN(length(bins), 1);
varV2 = NaN(length(bins), 1);
semV2 = NaN(length(bins), 1);
for i = 2 : length(bins)
    selectedV2 = Vart2(Vart1 <= bins(i) & Vart1 >= bins(i-1));
    meanV2(i-1) = nanmean(selectedV2);
    varV2(i-1) = nanvar(selectedV2);
    n = length(selectedV2);
    semV2(i-1) = nanstd(selectedV2)/sqrt(n);
    
    %         if ~isempty(selectedV2)
    %         [a, b] = ksdensity(selectedV2);
    %         plot(b, a, 'Color', [i/length(bins) i/length(bins) i/length(bins)], 'Linewidth', 2)
    %         hold on
    %         end
end

%***
fig = figure;
title(dt)
%scatter(Vart1(:), Vart2(:), 12)
hold on
plot(bins + w/2, meanV2, 'Linewidth', 1)
errorbar(bins + w/2, meanV2, semV2, '.', 'Color', [0.2 0.2 0.2])
%plot(bins + w/2, varV2, 'Linewidth', 1)
grid on
%xlim([prctile(Vart1(:),0.5) prctile(Vart1(:),99.5)])
%xlim([-0.2 0.2])
%ylim([-yl yl])
%xlabel(labelx)
%ylabel(labely)
ax = gca;
ax.FontSize = 14;


%% --- memory : multiple bouts ---


opt = 1
fig = figure;

for dt = 1
    
    if opt == 1
        Vart1 = dLu(:, 1:end-1);
        labelx = 'dLu';
        Vart2 = dX(:, 2:end).*(dX(:, 1:end-1)./abs(dX(:, 1:end-1)));
        labely = 'dX';
        yl=0.4;
    elseif opt == 2
        Vart1 = dLu(:, 1:end-dt);
        labelx = 'dLu';
        Vart2 = dLu(:, dt+1:end);
        labely = 'dLu';
        yl=0.1;
    elseif opt == 3
        Vart1 = dLu(:, 1:end-dt);
        Vart2 = dX(:, dt+1:end);
        labelx = 'dLu';
        labely = 'dXwrapped';
        yl=0.4;
    end
    
    b = max( abs(min(Vart1(:))), max(Vart1(:)));
    bins = linspace(-b, b, 40);
    w = mean(diff(bins));
    
    meanV2 = NaN(length(bins), 1);
    varV2 = NaN(length(bins), 1);
    semV2 = NaN(length(bins), 1);
    for i = 2 : length(bins)
        selectedV2 = Vart2(Vart1 <= bins(i) & Vart1 >= bins(i-1));
        meanV2(i-1) = nanmean(selectedV2);
        varV2(i-1) = nanvar(selectedV2);
        n = length(selectedV2);
        semV2(i-1) = nanstd(selectedV2)/sqrt(n);
        
%         if ~isempty(selectedV2)
%         [a, b] = ksdensity(selectedV2);
%         plot(b, a, 'Color', [i/length(bins) i/length(bins) i/length(bins)], 'Linewidth', 2)
%         hold on
%         end
    end
    
    %***
    %subplot(1,1,dt)
    title(dt)
    %scatter(Vart1(:), Vart2(:), 12)
    hold on
    plot(bins + w/2, meanV2, 'Linewidth', 1)
    errorbar(bins + w/2, meanV2, semV2, '.', 'Color', [0.2 0.2 0.2])
    %plot(bins + w/2, varV2, 'Linewidth', 1)
    grid on
    %xlim([prctile(Vart1(:),0.5) prctile(Vart1(:),99.5)])
    %xlim([-0.2 0.2])
    %ylim([-yl yl])
    xlabel(labelx)
    ylabel(labely)
    ax = gca;
    ax.FontSize = 14;
    
end
