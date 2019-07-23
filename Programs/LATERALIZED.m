%% ANALYSIS OF STEREOVISUAL PHOTOTAXIS DATA

%% ::: Load data :::
Lat.Load

%% Get corresponding date/fish from fishID (fishN)
figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Figures201808/';
path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
name = 'lateralized_info.mat';
load([path name], 'Linfo')
dat = Linfo{1};
fish = Linfo{2};
correspondance = Linfo{3};

fish_u_r_looking_for = 30;
correspmat = (correspondance == fish_u_r_looking_for);
correspfish = find(correspmat);
corresprow = find(sum(correspmat,2));
dat(corresprow)
fish(correspfish)
%% ::: basic stats :::
% some stats on sequences
    % number of fish
seqPerFish = histc(FishID, unique(FishID));
NFish = length(seqPerFish);
meanNofSeqPerFish = mean(seqPerFish);

NofSeqs = size(XLat, 1);
boutsPerSeq = size(XLat, 2) - sum(isnan(XLat),2);
maxnseq = max(seqPerFish);
medboutsperseq = median(boutsPerSeq);

fprintf('Fish N = %d .\n', NFish)
fprintf('Sequences: %d .\n', NofSeqs)
fprintf('Bouts/sequence median %f, mean %f .\n', median(boutsPerSeq), mean(boutsPerSeq))

%***
figure
subplot(1,2,1)
histogram(seqPerFish)
title('sequences/fish')
 
subplot(1,2,2)
histogram(boutsPerSeq)
title('bouts/sequence')
%% ::: Illumination profile (theoretical & on data) :::

Lat.illum_profile(Xfilt, L, R)


%% check coherence of variables

%***
figure
subplot(3,1,1)
plot(rad2deg(XLat), rad2deg(Xfilt), '.')
xlabel('XLat')
ylabel('Xfilt')

subplot(3,1,2)
err = XLat-Xfilt;
err(err==0)=[];
histogram(rad2deg(err))
title('error on angle detection in the experiment')

L2 = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),pi-abs(pi-mod(Xfilt,2*pi))); 
R2 = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),abs(pi-mod(Xfilt,2*pi)));
DIlr2 = L2-R2;
err_i = DIlr-DIlr2;
subplot(3,1,3)
histogram(err_i(err_i~=0))

% --- gaussian fit of angular error---
err(isnan(err))=[];
N = length(err);
binwidth = 3*iqr(err(:))/(N^(1/3));
bins = min(err(:)) : binwidth : max(err(:));
[yhist, edges]  = histcounts(err, bins);
xhist = edges(1:end-1) + binwidth/2;
yhist = yhist/(sum(yhist)*binwidth);
f = fit(xhist.',yhist.','gauss1');
sigma_error = f.c1;

%***
figure
plot(f, xhist, yhist)

%% ::: Basic visualization :::

edit Lat.basic_visu_perfish...(XLat, Xfilt, DIlr, FishN, visux, traj)

%% delete some sequences to get a uniform initial distribution
theta1 = wrapToPi(XLat(:,1));
[y, x] = hist(theta1, 12);
[~, maxidx] = max(smooth(y));
toomany = x(maxidx);

startatmax = find(theta1 > toomany-pi/3 & theta1 < toomany+pi/3);
randdel = ceil(rand(round(length(startatmax)/1.5),1)*length(startatmax-1));

Xranddel = XLat;
Xranddel(randdel,:) = [];

histogram(theta1)
hold on
histogram(wrapToPi(Xranddel(:,1)))
figure
polarhistogram(Xranddel(:,1))

%% ::: X distribution :::
Nbins=18;
colourID = 3;
[fig] = polar_distribution(XLat+pi/2, Nbins, FishID, colourID);

[figsub] = polar_distribution(XLat(:,2:17)+pi/2, Nbins, FishID, colourID);
[figsub] = polar_distribution(xfish+pi/2, Nbins, 1:20, colourID);


[figsub2] = polar_distribution(XLat(:,2:boutsPerSeq)+pi/2, Nbins, FishID, colourID);


% --- Some variants ---
%***
figure
subplot(2,2,1)
polarhistogram(Xfilt(:))
title('Biased - distribution of X at each bout')

subplot(2,2,2)
subX = Xfilt(:,3:20);
polarhistogram(subX(:))
title('bout# 3-20')

subplot(2,2,3)
polarhistogram(X(:))
title('time < 60 s')

subplot(2,2,4)
histogram(boutsPerSeq, 100)
title('bouts/sequence')
xlim([0 100])

%% ::: dX distribution :::
% Don't get wturn pturn from here : doesn't make sense
dX = dXl;
[wturn, wfor, pturn] = Spont.dXdistribution(dX, FishN, unique(FishN))


%% Plot various variables against eachother

v=1;
if v == 1
    %Vart1 = wrapToPi(XLat(:, 1 :end-1));
    Vart1 = DIlr(:, 1:end-1);
    Vart2 = dXl;
    b=20;
elseif v ==2
    Vart1 = DIlr(:, 1:end-1);
    Vart2 = Dist(:, 1:end)/pxmm;
    b = 30;
elseif v ==3
    Vart1 = DIlr(:, 1:end-1);
    Vart2 = Advancement(:, 1:end)/pxmm;   
    b=15;
elseif v ==4
    Vart2 = (dXl(:, 1:end-1).*dXl(:, 2:end))./(abs(dXl(:, 1:end-1)).*abs(dXl(:, 2:end)));
    Vart1 = IBI(:, 1:end-1);
    b=12;
elseif v ==5
    Vart1 = DIlr(:, 1:end-1);
    Vart2 = T(:, 1:end)/pxmm;
    b=15;
end

% --- bins with equal number of elements
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, b, b+3);
mV2 = nanmean(v2bin,2);
mV2sq = nanmean(v2bin.^2,2);

%***
figure;

subplot(2,1,1)
errorbar(binvals, mV2, std(v2bin,1,2)/sqrt(elts_per_bin),...
    'Linewidth', 2, 'DisplayName', 'Data')
legend
ax = gca;
ax.FontSize = 14;

subplot(2,1,2)
errorbar(binvals, mV2sq, std(v2bin.^2,1,2)/sqrt(elts_per_bin-1), 'Linewidth', 2)
ax = gca;
ax.FontSize = 14;

%% R in time / per fish
    edit Lat.R
    
%% AUTOCORRELATION : reinforcement/inhibition
    edit Lat.Autocorrelation

%% DIFFUSION
    edit Lat.Diffusion
    
%% bout bias ::: <dX> all data and double gaussian fit
    edit Lat.Drift

%% --
% --- regular bins ---
Nbins = 10;
Vart1 = IBI(:, 1:end-1);
Vart2 = (dXl(:, 1:end-1).*dXl(:, 2:end))./(abs(dXl(:, 1:end-1)).*abs(dXl(:, 2:end)));
[binvals, v2mean, v2std, v2eltspb] = regular_bins(Vart1, Vart2, Nbins);

%***
figure
errorbar(binvals, v2mean, v2std./sqrt(v2eltspb), 'Linewidth', 2)
title('regularly interspaced bins')

%% model : stationary distribution

interval = 0 : 0.05 : pi;
A = 0.11;
Dcoeff = 0.3;

Pdistr = exp(-A/Dcoeff.*interval.^2);
norm = trapz(interval, Pdistr);

%***
figure
subplot(2,1,1)
histogram(abs(wrapToPi(XLat(:,1:medboutsperseq)+pi/2)),interval, 'Normalization', 'pdf');
hold on
plot( interval, 0.7*Pdistr/norm)
xticks([0,pi/2, pi, 3*pi/2, 2*pi])
xticklabels({'0','\pi/2', '\pi'})
legend('exp(-\theta^2*A/D)')

subplot(2,1,2)
theta = wrapToPi(XLat(:,1:medboutsperseq)+pi/2);
polarhistogram(theta(:), 'Normalization', 'pdf');
hold on
polarplot([-interval NaN interval], [0.35*Pdistr/norm NaN 0.35*Pdistr/norm])

%% trajectory visualization
xCoord0 = (xCoord-xCoord(:,1))/pxmm;
yCoord0 = (yCoord-yCoord(:,1))/pxmm;
rho = sqrt( (xCoord0).^2 + (yCoord0).^2 );

x = xCoord0.*cos(-wrapToPi(XLat(:,1))+Xlab(:,1)) - yCoord0.*sin(-wrapToPi(XLat(:,1)) + Xlab(:,1));
y = xCoord0.*sin(-wrapToPi(XLat(:,1))+Xlab(:,1)) + yCoord0.*cos(-wrapToPi(XLat(:,1)) + Xlab(:,1));
plot(x',y','.-')
title('source to the top')


%***
figure
adv = Advancement;
polarplot(XLat, [1:size(XLat,2)])
