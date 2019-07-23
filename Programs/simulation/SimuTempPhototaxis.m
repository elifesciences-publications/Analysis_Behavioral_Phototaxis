
%% ORIENTATION
%%
% load data
path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
file = 'temporal_exps';
load([path file], 'E')

for i = 1 : size(E, 2)
    if E(i).ExpType == 'exp60'
        e6_x = E(i).AngleSource;
        e6_lum = E(i).Lum;
        e6_f = E(i).FishN;
        e6_r = E(i).R;
        e6_tb = E(i).TimeBout;
    elseif E(i).ExpType == 'exp30'        
        e3_x = E(i).AngleSource;
        e3_lum = E(i).Lum;
        e3_f = E(i).FishN;
        e3_r = E(i).R;
        e3_tb = E(i).TimeBout;
    elseif E(i).ExpType == 'sin60'
        s6_x = E(i).AngleSource;
        s6_lum = E(i).Lum;
        s6_f = E(i).FishN;
        s6_r = E(i).R;
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
lum_exp6 = luminosity_exponentielle(0.6);
lum_exp3 = luminosity_exponentielle(0.3);
lum_sin6 = luminosity_sinus(0.6);

[dx_histogram, x_histogram] = hist(alldx, 100);
f = fit(x_histogram.',dx_histogram.','gauss2'); 

%% check proba functions
dll = -2:0.01 : 1;
[pt, wt] = ProbaTurnFLum(dll);

%***
figure;
plot(dll, pt, dll, wt);
legend('proba turn', 'w turn')

%% choose experiment 
lum = lum_exp3;
x = e3_x;
f = e3_f;
boutsperseq = size(x,2)-sum(isnan(x),2);
medboutsperseq = median(boutsperseq);

%% parameters

%rng('shuffle')
Nexp = 10000; % number of experiments 
Nsteps = 50; % number of time steps per experiment

[theta_complete, luminosity_complete, lambda] = random_walk_temporal_bias(Nexp, Nsteps, lum);


alpha = theta_complete(:);
[mp, rho_p, mu_p] = circ_moment(alpha);

Nbins = 18;
angle = ( (1:Nbins) -0.5) * 2*pi/Nbins;
PDFangle = hist(mod(theta_complete(:),2*pi),Nbins)/sum(hist(mod(theta_complete(:),2*pi),Nbins));
PDFangleExp =  hist(mod(x(:),2*pi),Nbins)/sum(hist(mod(x(:),2*pi),Nbins));

% ***
fig = figure;

polarplot([angle,angle(1)],[PDFangle,PDFangle(1)], 'k', 'Linewidth', 2);
hold on
polarplot([angle,angle(1)],[PDFangleExp,PDFangleExp(1)], 'Linewidth', 1);
set(gca,'ThetaZeroLocation','top',...
    'ThetaDir','clockwise')
thetaticks( 0 :90: 360 )
thetaticklabels({'0' '\pi/2' '\pi' '-\pi/2'})
title({'sinusoidal 60%', 'X at each bout (10000 sequences, 30 bouts), X \in [-\pi \pi]'})
legend('simu, R = 0.16', 'experiment')
ax = gca;
ax.RAxisLocation = 45;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
ax.FontWeight = 'normal';

% calculate mean angle and moment of the distribution.
display(['mean angle in degree  :  ',num2str(mu_p*180/pi)]);
display(['first moment  :  ',num2str(rho_p)]);
[p,z] = circ_rtest(theta_complete(:)/180*pi);
disp(['p = ' num2str(p) ' ; z = ' num2str(z)])

disp('R proj')
circ_r(theta_complete(:))*cos(circ_mean(theta_complete(:)-pi)) 

% ***
med_bouts_per_seq = median(size(x,2)-sum(isnan(x),2));
xData = x(:,1:medboutsperseq);
xSimu = theta_complete;

colourID = 4;
[fig] = polar_distribution_simu_exp(xData+pi, xSimu+pi, Nbins, f, colourID)

toc

% --- ROI ---
virtual_roi_radius = 41 ;

% calculate cartesian coordinates at each time point
theta = theta_complete;
a = NaN(Nexp, Nsteps);
b = NaN(Nexp, Nsteps);
%a(:,1) = 0;
%b(:,1) = 0;
rho_ini = 20+1.3*randn(Nexp, 1);
theta_ini = deg2rad(randi(360, Nexp, 1));
[a(:,1),b(:,1)]= pol2cart(theta_ini, rho_ini);
%ind = randsample(length(a_ini), Nexp, true);
%a(:,1) = a_ini(ind)-41;
%b(:,1) = b_ini(ind)-41;

for i = 1 : Nsteps-1
    [da, db] = pol2cart(theta(:,i), lambda(:,i));
    a(:,i+1) = a(:,i) + da;
    b(:,i+1) = b(:,i) + db;
end

% finalize
within_roi = (sqrt(a.^2 + b.^2) < virtual_roi_radius);
theta(within_roi == 0) = NaN;
meanThovertimeSimu = circ_mean(theta-pi);
RcircSimu = circ_r(theta-pi);
RprojSimu = RcircSimu.*cos(meanThovertimeSimu);

boutsperseqSimu = size(theta,2)-sum(isnan(theta),2);
mednboutsSimu = median(boutsperseqSimu);

RcircData = circ_r(x-pi);
meanThovertimeData = circ_mean(x-pi);
RprojData = RcircData.*cos(meanThovertimeData);

%***
figure
plot(RprojSimu, 'LineWidth', 2)
hold on
plot(RprojData, 'LineWidth', 2)
xlim([0 30])
legend('simu', 'data')
xlabel('bout #')
ylabel('R')
ax=gca;
ax.FontSize = 16;
ax.LineWidth = 1.5;

colourID = 4;
%***
[fig] = polar_distribution_simu_exp(xData(:,1:medboutsperseq)+pi, theta(:,1:mednboutsSimu)+pi, Nbins, f, colourID)

%% trajectories
%***
seqmax = Nexp
aa = a(1:seqmax,:);
bb = b(1:seqmax,:);
th = theta(1:seqmax,:);
ll = lambda(1:seqmax,:);
wr = within_roi(1:seqmax,:);
aa(wr == 0) = NaN;
bb(wr == 0) = NaN;
for i = 1 : size(aa,1)
    seqend = find(isnan(aa(i,:)),1);
    if ~isempty(seqend)
        aa(i, seqend:end) = nan;
        bb(i, seqend:end) = nan;
        th(i, seqend:end) = nan;
    end
end

%plot(aa', bb', '.-')
hold on
%plot(aa(:,1)', bb(:,1)', 'sq')

ll(wr == 0) = NaN;

%***
%figure
xCoord0 = (aa-aa(:,1));
yCoord0 = (bb-bb(:,1));
rho = sqrt( (xCoord0).^2 + (yCoord0).^2 );

arot = xCoord0.*cos(-pi/2) - yCoord0.*sin(-pi/2);
brot = xCoord0.*sin(-pi/2) + yCoord0.*cos(-pi/2);
%plot(arot',brot','.-')
%title('source to the top')

%% Compare R data and simu stationary
xsub = x(:,2:medboutsperseq);
RcircData = circ_r(xsub(:)-pi);
meanTh = circ_mean(xsub(:)-pi);
RprojData = RcircData.*cos(meanTh)

medboutsperseqsimu = median(size(theta,2) - sum(isnan(theta),2));
thetasub = theta(:,2:medboutsperseqsimu);
meanThsimu = circ_mean(thetasub(:)-pi);
RcircSimu = circ_r(thetasub(:)-pi);
RprojSimu = RcircSimu.*cos(meanThsimu)
%% time evolution
b = 1 : 10 : Nsteps;
for i = 2 : Nsteps
    polarhistogram(theta_complete(:,i), 50)
    drawnow
    text(30, 15, num2str(i))
    pause(0.2)
end
figure
plot(circ_r(theta_complete, [], [], 1))

%%
dtheta = diff(theta_complete, 1, 2);
figure
subplot(1,2,1)
histogram(dtheta)
subplot(1,2,2)
plot(xcorr(dtheta(1,:)))

%% 
lum = luminosity_complete;
dlum = diff(luminosity_complete,1,2);
dx = diff(theta_complete,1,2);
v1 = dlum(:, 1:end-1)./( (lum(:,2:end-1)+lum(:,1:end-2))./2 );
v2 = dx(:, 2:end);
