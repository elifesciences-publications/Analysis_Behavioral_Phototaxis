%% Simulation sterovisual phototaxis

% --- Fetch data ---
if ~exist('XLat','var')
    path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
    name = 'lateralized_exps.mat';
    load([path name], 'El')
    thetaData = El.AngleSource;
    thetaDatafilt = El.AngleSourceFiltered;
    f = El.FishN;
    medboutsperseq = median(size(thetaData,2) - sum(isnan(thetaData),2));
    a_ini = El.xCoord(:,1)/11.5;
    b_ini = El.yCoord(:,1)/11.5;
    a_ini(isnan(a_ini))=[];
    b_ini(isnan(a_ini))=[];
    clear El
end

% === PARAMETERS TO TUNE ===
Nexp=500000; % number of independent experiments
Nsteps=50;  % number of time steps

w1=wturn;   % std of turning bouts amplitude distribution 
w2=wfor;    % std of forward bouts amplitude distribution 
p=pturn;    % probability of turning bout
pflip=0.19;  % neutral probability of flipping orientation
%boutn = 1:Nsteps;
%bias_slope = 20*(0.22-0.1*(log(-boutn/4))); % ratio of bias to contrast
%bias_slope(bias_slope<0) = 0;
bias_slope = 6;% 6 if decremental, 2.3 if constant

% --- initiate values ---
lrst=(rand(Nexp,1)>0.5)*2-1;    % initial left/right state
dx=[];                          % change in orientation
theta=-(rand(Nexp, 1)-0.5)*2*pi;    % uniform initial orientation

ind = randsample(length(thetaData(:,1)), Nexp, true);
theta = thetaData(ind,1)+pi/2;             % sample from data
theta = rand(Nexp,1)*2*pi;

qq=[];                          % flipping probability
lr=[];                          % left right state
lambda=[];
TS = zeros(Nexp, Nsteps);

for t=1:Nsteps
    xlast=wrapToPi(theta(:,end));
    q = pflip*(1-0.5*bias_slope/(t^.8)*lrst.*stereo_bias_mat(xlast));
    %q = pflip*(1-0.5*bias_slope*lrst.*stereo_bias_mat(xlast));
    q(q<0)=0;
    q(q>1)=1;
    lrst=lrst.*((rand(Nexp,1)>q)*2-1);
    turn=rand(Nexp,1)<p;
    dxt=zeros(Nexp,1);
    indturn=find(turn==1);
    indfor=find(turn==0);
    dxt(indturn)=lrst(indturn).*abs(w1*randn(length(indturn),1));
    dxt(indfor)=w2*randn(length(indfor),1);
    dx=[dx,dxt];
    theta=[theta,theta(:,end)+dxt];
    qq=[qq, q];
    lr=[lr, lrst];
    l=NaN(Nexp,1);
    l(indturn) = gamrnd(3.5, 1, length(indturn),1);
    l(indfor) = gamrnd(2.8, 1, length(indfor),1);
    lambda = [lambda, l];
    TS(:,t) = turn;
end

%% --- ROI ---
virtual_roi_radius = 41 ;

% calculate cartesian coordinates at each time point
thetaSimuROI = theta(:,1:end-1);
a = NaN(Nexp, Nsteps);
b = NaN(Nexp, Nsteps);
%a(:,1) = 0;
%b(:,1) = 0;
rho_ini = 20+1.3*randn(Nexp, 1);
theta_ini = deg2rad(randi(360, Nexp, 1));
[a(:,1),b(:,1)]= pol2cart(theta_ini, rho_ini);
ind = randsample(length(a_ini), Nexp, true);
a(:,1) = a_ini(ind)-41;
b(:,1) = b_ini(ind)-41;

for i = 1 : Nsteps-1
    [da, db] = pol2cart(thetaSimuROI(:,i), lambda(:,i));
    a(:,i+1) = a(:,i) + da;
    b(:,i+1) = b(:,i) + db;
end

% finalize
within_roi = (sqrt(a.^2 + b.^2) < virtual_roi_radius);
thetaSimuROI(within_roi == 0) = NaN;
aROI = a;
bROI = b;
aROI(within_roi == 0) = NaN;
bROI(within_roi == 0) = NaN;

% delete all poits after an exist
thetanans = isnan(thetaSimuROI);
[~, firstnancol] = max(thetanans, [], 2 );
for i = 1 : size(aROI,1)
    if ~isempty(seqend)
        aROI(i, firstnancol(i):end) = nan;
        bROI(i, firstnancol(i):end) = nan;
        thetaSimuROI(i, firstnancol(i):end) = nan;
    end
end

% get scattering R
meanThovertime = circ_mean(thetaSimuROI);
Rcirc = circ_r(thetaSimuROI);
Rproj = Rcirc.*cos(meanThovertime);

boutsperseqsimu = size(thetaSimuROI,2)-sum(isnan(thetaSimuROI),2);
mednboutssimu = median(boutsperseqsimu);

% --- get corresponding contrast ---
percPm =1;
lum_lin = luminosity_linear(percPm);
L = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),pi-abs(pi-mod(theta-pi/2,2*pi))); 
R = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),abs(pi-mod(theta-pi/2,2*pi))); % experiment shifted by pi/2
Csimu = L-R;
Csimu = Csimu/max(abs(Csimu(:)));

% --- data R ---
meanThovertime = circ_mean(thetaData+pi/2);
RcircData = circ_r(thetaData+pi/2);
RprojData = RcircData.*cos(meanThovertime);

boutsperseq = size(thetaData,2)-sum(isnan(thetaData),2);
mednbouts = median(boutsperseq);

%... per fish
% data
nbouts = 50;
id = unique(f);
RcircDataIndiv = nan(length(id), nbouts);
meanThovertimeDataIndiv = nan(length(id), nbouts);
for i = 1 : length(id)
    indiv = find(f == i);
    thetaDataSample = thetaData(indiv,1:nbouts);
    RcircDataIndiv(i,:) = circ_r(thetaDataSample+pi/2);
    meanThovertimeDataIndiv(i,:) = circ_mean(thetaDataSample+pi/2);
end
RprojDataIndiv = RcircDataIndiv.*cos(meanThovertimeDataIndiv);

%***
figure
plot(Rproj)
hold on
plot(RprojData, 'k', 'Linewidth', 1.5)
errorbar([1:nbouts], mean(RprojDataIndiv,1), std(RprojDataIndiv,1)/sqrt(length(id)))
xlim([0 50])

% display a sample
%***
figure
seqmax = 1000;
plot(a(1:seqmax,:)', b(1:seqmax,:)', '.-')
hold on
plot(a(:,1)', b(:,1)', 'sq')

% trajectories from 0
xCoord0 = (aROI-aROI(:,1));
yCoord0 = (bROI-bROI(:,1));
rho = sqrt( (xCoord0).^2 + (yCoord0).^2 );

arot = xCoord0.*cos(pi/2) - yCoord0.*sin(pi/2);
brot = xCoord0.*sin(pi/2) + yCoord0.*cos(pi/2);

%***
figure
plot(arot(1:seqmax,:)', brot(1:seqmax,:)', '.-')
hold on
plot(arot(:,1)', brot(:,1)', 'sq')
title('source to the top')

%% save trajectory

% save trajectory
traj.exptype = 'no_stim';
traj.x = arot;
traj.y = brot;

save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/AnalysisOutput/trajectories/' traj.exptype '.mat'], 'traj')


%%
colourID = 4;
Nbins = 18;
fig = polar_distribution_simu_exp(thetaData(:,1:medboutsperseq)+pi/2, theta, Nbins, f, colourID);
title('without escaping ROI')

fig2 = polar_distribution_simu_exp(thetaData(:,1:medboutsperseq)+pi/2, thetaSimuROI(:,1:mednboutssimu), Nbins, f, colourID);
title('with escaping ROI')

Nbins = 18;
angle = ((1:Nbins) - 0.5) * 2*pi/Nbins;
xpdf_data = hist(mod(thetaData(:),2*pi),Nbins)/sum(hist(mod(thetaData(:),2*pi),Nbins));
XLsub = thetaData(:,1:medboutsperseq);
xpdf_datasub = hist(mod(XLsub(:),2*pi),Nbins)/sum(hist(mod(XLsub(:),2*pi),Nbins));
xpdf_simu = hist(mod(theta(:)-pi/2,2*pi),Nbins)/sum(hist(mod(theta(:)-pi/2,2*pi),Nbins));

%***
figure
polarplot(angle, xpdf_datasub, 'Linewidth', 1.5, 'DisplayName', 'subdata')
hold on
polarplot(angle,xpdf_data, 'Linewidth', 1.5, 'DisplayName', 'data')
polarplot(angle,xpdf_simu, 'Linewidth', 1.5, 'DisplayName', 'simu')
legend

% --- extract the drift ---
xwrap2pi=wrapToPi(theta(:,1:end-1));
[binvals,elts_per_bin, v2bin, ~, ~ ] = ...
    BinsWithEqualNbofElements(Csimu(:,1:end-1), dx, 10, 12 );
Drift = mean(v2bin,2);
pol = polyfit(binvals,Drift,1);
Drift_fit= polyval(pol,binvals);

[binvals_sub,elts_per_bin_sub, v2bin_sub, ~, ~ ] = ...
    BinsWithEqualNbofElements(Csimu(:,1:25), dx(:,1:25), 10, 12 );
Drift_sub = mean(v2bin_sub,2);
pol_sub = polyfit(binvals_sub,Drift_sub,1);
Drift_sub_fit= polyval(pol_sub,binvals_sub);

figure
subplot(1,2,2)
plot(binvals,Drift);
hold on;
plot(binvals, Drift_fit)
%DD=[DD, p(1)]
%end
pflip*bias_slope*sqrt(2/pi)*w1*p/pol(1)
%ylim([-0.1,0.1])
hold on;
plot([-1 1], [0 0])
plot([0 0] , [-.1 .1])

subplot(1,2,1)
plot(binvals_sub,Drift_sub);
hold on;
plot(binvals_sub, Drift_sub_fit)

%% plot histogram
tinit=5;
tend=11;
xpart=abs(wrapToPi(theta(:,tinit:tend)));
histogram(xpart(:),15, 'Normalization', 'pdf');

%%

Nexp_real=300; % number of trajectories in actual experiments;
for p=1:1000;
    xpart=theta(fix(Nexp*abs(rand(Nexp_real,1)))+1,tinit:tend);
    hval(p,:)=histcounts(abs(wrapToPi(xpart(:))),binedges, 'Normalization', 'pdf');
    %hval(p,:)=hpart.Values;
end
hstd=std(hval,1);
figure;
boundedline([1:15]/15*pi,hmean_val, hstd);

%% extract the drift
dxt = dx;
dxt(TS==0) = NaN;
dxt(dxt<0) = -1; dxt(dxt>0) = 1;
ac1 = dxt(:,1:end-1).*dxt(:,2:end);
[binvals,elts_per_bin, v2bin, ~, ~ ] = ...
    BinsWithEqualNbofElements(Csimu(:,1:end-2), ac1, 10, 12 );
mac1 = mean(v2bin,2);
figure;
plot(binvals, mac1)

%% get effective diffusion coefficient

dxbias=pol(1)*wrapToPi(theta(:,1:end-1));
plot(wrapToPi(theta(:,1:end-1)), dxbias);
xdiff=cumsum(dx+dxbias);
MSD=[];
for n=1:max(Nsteps-1,100);
    D2=diff(xdiff(1:n:end,:),1).^2;
    MSD(n)=mean(D2(:));
end
% figure;
% plot([0:max(Nsteps-1,100)],[0 MSD])
% fitrange=[5:max(Nsteps-1,100)];
% Deff_fit= polyfit(fitrange,MSD(fitrange),1);
% Deff=Deff_fit(1);
%% Compare with O.U. process
figure;
histogram(abs(wrapToPi(theta)),30, 'Normalization', 'pdf')
k=pol(1)/Deff/2;
xfit=[0:pi/100:pi];
norm=1/2*sqrt(pi/k)*erf(sqrt(k)*pi);
hfit=1/norm*exp(-k*xfit.^2);
hold on;
plot(xfit,hfit)

%xlim([0,pi/4])
%%
ths=thetaSimuROI(:,1:mednboutssimu);
thd=thetaData(:,1:medboutsperseq)+pi/2;

ths = ths(:);
ths(isnan(ths))=[];
thd = thd(:);
thd(isnan(thd))=[];

figure
histogram(wrapToPi(ths),36, 'Normalization', 'pdf')
hold on
histogram(wrapToPi(thd),36, 'Normalization', 'pdf')
ax = gca;
bwidth = ax.Children(1).BinWidth;
xhist = ax.Children(1).BinEdges(1:end-1) + bwidth ;

k = 0.06/0.3; % A/D
xfit=[-pi xhist pi];
hfit=exp(-k*(xfit).^2);
norm=sqrt(k/pi);
f = @(x) exp(-k*x.^2);
hfitnorm = hfit/integral(f, -pi, pi);
plot(xfit,hfitnorm, 'LineWidth', 1.5)
legend('simulation', 'data', 'analytical solution')

% circular
ths = ths(:);
ths(isnan(ths))=[];
thd = thd(:);
thd(isnan(thd))=[];

figure
polarhistogram(wrapToPi(ths),36, 'Normalization', 'pdf')
hold on
polarhistogram(wrapToPi(thd),36, 'Normalization', 'pdf')
polarplot(xfit,hfitnorm, 'LineWidth', 1.5)

%%
%dxcropped = dx(:,1:end-1);
dxcropped = dx;
dxcropped(TS==0)=NaN;
xwrap2pi=wrapToPi(theta);
Lat.a.auto_co_reinforcement(sign(dxcropped), Csimu, 20, 15)

%% check inter-bout distance

[binvals,elts_per_bin, v2bin, ~, ~ ] = ...
    BinsWithEqualNbofElements(dx, lambda, 20, 22 );
meanlambda = mean(v2bin,2);
varlambda = var(v2bin,1,2);
figure
plot(binvals, meanlambda);
hold on
plot(binvals, varlambda);

figure
histogram(lambda(TS==1), 'Normalization', 'pdf')
hold on
histogram(lambda(TS==0), 'Normalization', 'pdf')