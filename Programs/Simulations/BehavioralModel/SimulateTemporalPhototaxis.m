
%% ORIENTATION

%% Load experimental data (once)
path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
file = 'temporal_exps';
load([path file], 'E')

for i = 1 : size(E, 2)
    if E(i).ExpType == 'exp60'
        e6_x = E(i).AngleSource;
        e6_lum = E(i).Lum;
        e6_f = E(i).FishN;
        e6_r = E(i).R;
        e6_d = E(i).Dist;
        e6_tb = E(i).TimeBout;
    elseif E(i).ExpType == 'exp30'        
        e3_x = E(i).AngleSource;
        e3_lum = E(i).Lum;
        e3_f = E(i).FishN;
        e3_r = E(i).R;
        e3_d = E(i).Dist;
        e3_tb = E(i).TimeBout;
    elseif E(i).ExpType == 'sin60'
        s6_x = E(i).AngleSource;
        s6_lum = E(i).Lum;
        s6_f = E(i).FishN;
        s6_r = E(i).R;
        s6_d = E(i).Dist;
        s6_tb = E(i).TimeBout;
    end
end

% dlum
e3_dlum = diff( e3_lum, 1, 2 );
e6_dlum = diff( e6_lum, 1, 2 );
s6_dlum = diff( s6_lum, 1, 2 );

% dx
e3_dx =  diff( e3_x, 1, 2 );
e6_dx =  diff( e6_x, 1, 2 );
s6_dx =  diff( s6_x, 1, 2 );
alldx = [e3_dx(:) ; e6_dx(:) ; s6_dx(:)];

% lum
lum_profile_exp6 = luminosity_exponentielle(0.6);
lum_profile_exp3 = luminosity_exponentielle(0.3);
lum_profile_sin6 = luminosity_sinus(0.6);

alldx(isnan(alldx))=[];
[dx_histogram, x_histogram] = hist(alldx, 100);
f = fit(x_histogram.',dx_histogram.','gauss2'); 

% decrease in lum before first bout
lum_max = lum_profile_exp3(end,2)/0.3/1000;
ini_lum_change_s6 = s6_lum(:,1) - lum_max;
ini_lum_change_e6 = e6_lum(:,1) - lum_max;
ini_lum_change_e3 = e3_lum(:,1) - lum_max;
ini_rellum_change_s6 = (s6_lum(:,1) - lum_max)./(s6_lum(:,1) + lum_max);
ini_rellum_change_e6 = (e6_lum(:,1) - lum_max)./(e6_lum(:,1) + lum_max);
ini_rellum_change_e3 = (e3_lum(:,1) - lum_max)./(e3_lum(:,1) + lum_max);

%***
% figure
% 
% subplot(2,2,[1 3])
% hold on
% histogram(s6_lum, 30, 'Normalization', 'pdf')
% histogram(e6_lum, 30, 'Normalization', 'pdf')
% histogram(e3_lum, 30, 'Normalization', 'pdf')
% 
% subplot(2,2,2)
% histogram(ini_lum_change_s6, 30, 'Normalization', 'pdf')
% hold on
% histogram(ini_lum_change_e6, 30, 'Normalization', 'pdf')
% histogram(ini_lum_change_e3, 30, 'Normalization', 'pdf')
% 
% subplot(2,2,4)
% histogram(ini_rellum_change_s6, 30, 'Normalization', 'pdf', 'DisplayName', 's6')
% hold on
% histogram(ini_rellum_change_e6, 30, 'Normalization', 'pdf', 'DisplayName', 'e6')
% histogram(ini_rellum_change_e3, 30, 'Normalization', 'pdf', 'DisplayName', 'e3')
% legend
%% choose experiment 
unif = 1;

if unif == 1
    disp('unif 1 : sin60 chosen')
    lum_profile = lum_profile_sin6;
    lumData = s6_lum;
    inidlum = ini_lum_change_s6;
    inireldlum = ini_rellum_change_s6;
    thetaData = s6_x;
    f = s6_f;
    exptype = 'unif_1';
elseif unif == 2
    disp('unif 2 : exp60 chosen')
    lum_profile = lum_profile_exp6;
    lumData = e6_lum;
    inidlum = ini_lum_change_e6;
    inireldlum = ini_rellum_change_e6;
    thetaData = e6_x;
    f = e6_f;
    exptype = 'unif_2';
elseif unif == 3
    disp('unif 3 : exp30 chosen')
    lum_profile = lum_profile_exp3;
    inidlum = ini_lum_change_e3;
    inireldlum = ini_rellum_change_e3;
    lumData = e3_lum;
    thetaData = e3_x;
    f = e3_f;
    exptype = 'unif_3';
end

boutsperseq = size(thetaData,2)-sum(isnan(thetaData),2);
medboutsperseqdata = median(boutsperseq);
reldlumData = diff(lumData(:,1:end-1),1,2)./((lumData(:,1:end-2) + lumData(:,2:end-1))/2);

%% check proba functions (just for display)
dll = -2:0.01 : 2;
[pt, wt] = ProbaTurnFLum(dll);

%***
figure;
plot(dll, pt, dll, wt);
ylim([0 max(wt)])
legend('proba turn', 'w turn')

%% Simulation

rng('shuffle')
Nexp = 10000; % number of experiments 
Nsteps = 100; % number of time steps per experiment

%ind = randsample(length(thetaData(:,1)), Nexp, true);
%theta_ini = wrapTo2Pi(thetaData(ind,1)); 
theta_ini = rand(Nexp,1)*2*pi;

[thetaSimu, luminositySimu, lambdaSimu] = random_walk_temporal_bias(Nexp, Nsteps, lum_profile, theta_ini);

%% histograms of luminosity profiles
lum_max_chosen = 0.05;
ini_lum_change_simu = luminositySimu(:,1)/1000  - lum_max_chosen;
ini_rellum_change_simu = (luminositySimu(:,1)/1000  - lum_max_chosen)./(luminositySimu(:,1)/1000  + lum_max_chosen);
reldlumSimu = diff(luminositySimu(:,1:end-1),1,2)./(luminositySimu(:,1:end-2)+luminositySimu(:,2:end-1))/2;

%***
figure

%data
subplot(2,2,[1])
hold on
lumdata_nonan = lumData(:);
lumdata_nonan(isnan(lumdata_nonan)) = [];
histogram(lumdata_nonan, 30, 'Normalization', 'pdf')

subplot(2,2,2)
hold on
reldlumData_nn = reldlumData(:);
reldlumData_nn(isnan(reldlumData_nn))= [];
histogram(reldlumData_nn, 30, 'Normalization', 'pdf')

subplot(2,2,3)
inidlum_nonan = inidlum(:);
inidlum_nonan(isnan(inidlum_nonan)) = [];
histogram(inidlum, 30, 'Normalization', 'pdf')

subplot(2,2,4)
inireldlum_nn = inireldlum(:);
inireldlum_nn(isnan(inireldlum_nn))=[];
histogram(inireldlum_nn, 30, 'Normalization', 'pdf', 'DisplayName', 's6')

% simulation
subplot(2,2,1)
hold on
histogram(luminositySimu/1000, 30, 'Normalization', 'pdf')
title('all luminosity')

subplot(2,2,2)
hold on
histogram(reldlumSimu, 30, 'Normalization', 'pdf')
title('all rel change intensity')

subplot(2,2,3)
hold on
histogram(ini_lum_change_simu, 30, 'Normalization', 'pdf')
title('initial lum change')

subplot(2,2,4)
hold on
histogram(ini_rellum_change_simu, 30, 'Normalization', 'pdf', 'DisplayName', 'simu')
title('initial rel lum change')


%% Distributions

[mp, rho_p, mu_p] = circ_moment(thetaSimu(:));

Nbins = 18;
angle = ( (1:Nbins) -0.5) * 2*pi/Nbins;
PDFangle = hist(mod(thetaSimu(:),2*pi),Nbins)/sum(hist(mod(thetaSimu(:),2*pi),Nbins));
PDFangleExp =  hist(mod(thetaData(:),2*pi),Nbins)/sum(hist(mod(thetaData(:),2*pi),Nbins));

% ***
fig = figure;
polarplot([angle,angle(1)],[PDFangle,PDFangle(1)], 'k', 'Linewidth', 2);
hold on
polarplot([angle,angle(1)],[PDFangleExp,PDFangleExp(1)], 'Linewidth', 1);
set(gca,'ThetaZeroLocation','top',...
    'ThetaDir','clockwise')
thetaticks( 0 :90: 360 )
thetaticklabels({'0' '\pi/2' '\pi' '-\pi/2'})
ax = gca;
ax.RAxisLocation = 45;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'normal';

% calculate mean angle and moment of the distribution.
display(['mean angle in degree  :  ',num2str(mu_p*180/pi)]);
display(['first moment  :  ',num2str(rho_p)]);
[p,z] = circ_rtest(thetaSimu(:)/180*pi);
disp(['p = ' num2str(p) ' ; z = ' num2str(z)])

disp('R proj')
circ_r(thetaSimu(:))*cos(circ_mean(thetaSimu(:)-pi)) 

% ***
Nbins = 18;
colourID = 4;
[fig] = polar_distribution_simu_exp(thetaData(:,1:medboutsperseqdata)+pi, thetaSimu+pi, Nbins, f, colourID);
title('no escape')

%% --- ROI ---
virtual_roi_radius = 41 ;

% calculate cartesian coordinates at each time point
a = NaN(Nexp, Nsteps);
b = NaN(Nexp, Nsteps);
%a(:,1) = 0;
%b(:,1) = 0;
rho_ini = 20+1.3*randn(Nexp, 1); %20
alpha_ini = deg2rad(randi(360, Nexp, 1));
[a(:,1),b(:,1)]= pol2cart(alpha_ini, rho_ini);
%ind = randsample(length(a_ini), Nexp, true);
%a(:,1) = a_ini(ind)-41;
%b(:,1) = b_ini(ind)-41;

for i = 1 : Nsteps-1
    [da, db] = pol2cart(thetaSimu(:,i), lambdaSimu(:,i));
    a(:,i+1) = a(:,i) + da;
    b(:,i+1) = b(:,i) + db;
end

% finalize
thetaROI = thetaSimu;
within_roi = (sqrt(a.^2 + b.^2) < virtual_roi_radius);
thetaROI(within_roi == 0) = NaN;
thetanans = isnan(thetaROI);
[~, firstnancol] = max( thetanans, [], 2 );
for i = 1 : size(thetaROI,1)
    thetaROI(i, firstnancol(i):end) = nan;
    a(i, firstnancol(i):end) = nan;
    b(i, firstnancol(i):end) = nan;
end
meanThovertimeSimu = circ_mean(thetaROI-pi);
RcircSimu = circ_r(thetaROI-pi);
RprojSimu = RcircSimu.*cos(meanThovertimeSimu);

boutsperseqSimu = size(thetaROI,2)-sum(isnan(thetaROI),2);
medboutsperseqSimu = median(boutsperseqSimu);

% data
nbouts = 30;
id = unique(f);
RcircDataIndiv = nan(length(id), nbouts);
meanThovertimeDataIndiv = nan(length(id), nbouts);
for i = 1 : length(id)
    indiv = find(f == i);
    thetaDataSample = thetaData(indiv,1:nbouts);
    RcircDataIndiv(i,:) = circ_r(thetaDataSample-pi);
    meanThovertimeDataIndiv(i,:) = circ_mean(thetaDataSample-pi);
end
RprojDataIndiv = RcircDataIndiv.*cos(meanThovertimeDataIndiv);

RcircData = circ_r(thetaData-pi);
meanThovertimeData = circ_mean(thetaData-pi);
RprojData = RcircData.*cos(meanThovertimeData);

%***
figure
plot(RprojSimu, 'LineWidth', 2, 'DisplayName', ['simu ' exptype])
hold on
plot(RprojData, 'LineWidth', 2, 'DisplayName', ['data ' exptype])
errorbar([1:nbouts], mean(RprojDataIndiv,1), std(RprojDataIndiv,1)/sqrt(length(id)))
xlim([0 30])
xlabel('bout #')
ylabel('R')
ax=gca;
ax.FontSize = 16;
ax.LineWidth = 1.5;

colourID = 4;
Nbins = 18;

%***
[fig] = polar_distribution_simu_exp...
    (thetaData(:,1:medboutsperseqdata)+pi, thetaROI(:,1:medboutsperseqSimu)+pi, Nbins, f, colourID);
title(['data med : ' num2str(medboutsperseqdata) ', ' 'simu med : ' num2str(medboutsperseqSimu)])

%% various parameters
%var1 = wrapToPi(thetaSimu(:,1:medboutsperseqSimu)-pi);
%var1 = (luminositySimu(:,1:end-1));
var1 = diff(luminositySimu(:,1:end-1),1,2)./(luminositySimu(:,1:end-2)+luminositySimu(:,2:end-1))/2;
var2 = diff(thetaSimu,1,2);
var2 = var2(:,2:end).^2;
%var2 = lambdaSimu;

bins = 200;

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, bins, bins+10);

%***
f = figure;
%scatter(var1(:), var2(:))
hold on
errorbar(binvals, mean(v2bin, 2), std(v2bin,1,2)/sqrt(elts_per_bin), 'Linewidth', 1.5)

%% trajectories
%***
seqmax = 100;

xCoord0 = (a-a(:,1));
yCoord0 = (b-b(:,1));
rho = sqrt( (xCoord0).^2 + (yCoord0).^2 );

arot = xCoord0.*cos(-pi/2) - yCoord0.*sin(-pi/2);
brot = xCoord0.*sin(-pi/2) + yCoord0.*cos(-pi/2);

%***
figure
plot(arot(1:seqmax,:)',brot(1:seqmax,:)','.-')
title('source to the top')

% save trajectory
traj.exptype = exptype;
traj.x = arot;
traj.y = brot;

save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/AnalysisOutput/trajectories/' exptype '.mat'], 'traj')
%% Compare R data and simu stationary
start_bout = 2
thsub = thetaData(:,start_bout:medboutsperseqdata);
RcircData = circ_r(thsub(:)-pi);
meanTh = circ_mean(thsub(:)-pi);
RprojData = RcircData.*cos(meanTh)

thetasub = thetaROI(:,start_bout:medboutsperseqSimu);
meanThsimu = circ_mean(thetasub(:)-pi);
RcircSimu = circ_r(thetasub(:)-pi);
RprojSimu = RcircSimu.*cos(meanThsimu)
%% time evolution
b = 1 : 10 : Nsteps;
for i = 2 : Nsteps
    polarhistogram(thetaSimu(:,i), 50)
    drawnow
    text(30, 15, num2str(i))
    pause(0.2)
end
figure
plot(circ_r(thetaSimu, [], [], 1))

%%
dtheta = diff(thetaSimu, 1, 2);
figure
subplot(1,2,1)
histogram(dtheta)
subplot(1,2,2)
autocorr = xcorr(dtheta);
plot(mean(autocorr(:,size(dtheta)+1:end)))

%% 
lum_profile = luminositySimu;
dlum = diff(luminositySimu,1,2);
dx = diff(thetaSimu,1,2);
v1 = dlum(:, 1:end-1)./( (lum_profile(:,2:end-1)+lum_profile(:,1:end-2))./2 );
v2 = dx(:, 2:end);
