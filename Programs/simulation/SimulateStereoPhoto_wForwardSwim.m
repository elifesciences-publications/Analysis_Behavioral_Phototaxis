%% Simulation sterovisual phototaxis

% === PARAMETERS TO TUNE ===
Nexp=100000; % number of independent experiments
Nsteps=50;  % number of time steps

w1=wturn;   % std of turning bouts amplitude distribution 
w2=wfor;    % std of forward bouts amplitude distribution 
p=pturn;    % probability of turning bout
pflip=0.188;  % neutral probability of flipping orientation
%boutn = 1:Nsteps;
%bias_slope = 20*(0.22-0.1*(log(-boutn/4))); % ratio of bias to contrast
%bias_slope(bias_slope<0) = 0;
bias_slope = 2;

% --- generate two populations of inter-bout distances ---
sigmafwd = 0.6;
sigmaturn = 0.5;
muturn=1;
mufwd = 0.7;
possible_lambda_values = 0:0.01:40;

% --- initiate values ---
lrst=(rand(Nexp,1)>0.5)*2-1;    % initial left/right state
dx=[];                          % change in orientation
x=-(rand(Nexp, 1)-0.5)*2*pi;    % uniform initial orientation
qq=[];                          % flipping probability
lr=[];                          % left right state
lambda=[];
TS = zeros(Nexp, Nsteps);

for t=1:Nsteps
    xlast=wrapToPi(x(:,end));
    q=pflip*(1-0.5*bias_slope*lrst.*stereo_bias_mat(xlast));
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
    x=[x,x(:,end)+dxt];
    qq=[qq, q];
    lr=[lr, lrst];
    l=NaN(Nexp,1);
    l(indturn) = pick_from_lognorm_distrib(possible_lambda_values, sigmaturn, muturn, length(indturn))';
    l(indfor) = pick_from_lognorm_distrib(possible_lambda_values, sigmafwd, mufwd, length(indfor))';
    lambda = [lambda, l];
    TS(:,t) = turn;
end

% --- ROI ---
% calculate cartesian coordinates at each time point
theta = x(:,1:end-1);
a = NaN(Nexp, Nsteps);
b = NaN(Nexp, Nsteps);
for i = 1 : Nsteps
    if i == 1
        a(:,i) = 0;
        b(:,i) = 0;
    else
        [da, db] = pol2cart(theta(:,i-1), lambda(:,i));
        a(:,i) = a(:,i-1) + da;
        b(:,i) = b(:,i-1) + db;
    end
end

% finalize
virtual_roi_radius = 42 ;
within_roi = (sqrt(a.^2 + b.^2) < virtual_roi_radius);
theta(within_roi == 0) = NaN;
meanThovertime = circ_mean(theta);
Rcirc = circ_r(theta);
Rproj = Rcirc.*cos(meanThovertime);
plot(Rproj)

% --- get corresponding contrast ---
percPm =1;
lum_lin = luminosity_linear(percPm);
L = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),pi-abs(pi-mod(x-pi/2,2*pi))); 
R = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),abs(pi-mod(x-pi/2,2*pi))); % experiment shifted by pi/2
Csimu = L-R;
Csimu = Csimu/max(abs(Csimu(:)));

% --- polar distribution of x ---
if ~exist('XLat','var')
    path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
    name = 'lateralized_exps.mat';
    load([path name], 'El')
    XLat = El.AngleSource;
    Xfilt = El.AngleSourceFiltered;
    FishID = El.FishN;
    medboutsperseq = median(size(XLat,2) - sum(isnan(XLat),2));
    clear El
end

colourID = 4;
fig = polar_distribution_simu_exp(XLat(:,1:medboutsperseq)+pi/2, x, Nbins, FishID, colourID);
title('without escaping ROI')

fig2 = polar_distribution_simu_exp(XLat(:,1:medboutsperseq)+pi/2, theta(:,1:30), Nbins, FishID, colourID);
title('with escaping ROI')

Nbins = 32;
angle = ((1:Nbins) - 0.5) * 2*pi/Nbins;
xpdf_data = hist(mod(XLat(:),2*pi),Nbins)/sum(hist(mod(XLat(:),2*pi),Nbins));
XLsub = XLat(:,1:medboutsperseq);
xpdf_datasub = hist(mod(XLsub(:),2*pi),Nbins)/sum(hist(mod(XLsub(:),2*pi),Nbins));
xpdf_simu = hist(mod(x(:)-pi/2,2*pi),Nbins)/sum(hist(mod(x(:)-pi/2,2*pi),Nbins));

%***
figure
polarplot(angle, xpdf_datasub, 'Linewidth', 1.5, 'DisplayName', 'subdata')
hold on
polarplot(angle,xpdf_data, 'Linewidth', 1.5, 'DisplayName', 'data')
polarplot(angle,xpdf_simu, 'Linewidth', 1.5, 'DisplayName', 'simu')
legend

% --- extract the drift ---
xwrap2pi=wrapToPi(x(:,1:end-1));
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
xpart=abs(wrapToPi(x(:,tinit:tend)));
histogram(xpart(:),15, 'Normalization', 'pdf');

%%

Nexp_real=300; % number of trajectories in actual experiments;
for p=1:1000;
    xpart=x(fix(Nexp*abs(rand(Nexp_real,1)))+1,tinit:tend);
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
ac1 = dxt(:,1:end-1).*dxt(:,2:end)
[binvals,elts_per_bin, v2bin, ~, ~ ] = ...
    BinsWithEqualNbofElements(Csimu(:,1:end-2), ac1, 10, 12 );
mac1 = mean(v2bin,2);
figure;
plot(binvals, mac1)

%% get effective diffusion coefficient

dxbias=pol(1)*wrapToPi(x(:,1:end-1));
plot(wrapToPi(x(:,1:end-1)), dxbias);
xdiff=cumsum(dx+dxbias);
MSD=[];
for n=1:max(Nsteps-1,100);
    D2=diff(xdiff(1:n:end,:),1).^2;
    MSD(n)=mean(D2(:));
end
figure;
plot([0:max(Nsteps-1,100)],[0 MSD])
fitrange=[5:max(Nsteps-1,100)];
Deff_fit= polyfit(fitrange,MSD(fitrange),1);
Deff=Deff_fit(1);
%% Compare with O.U. process
figure;
histogram(abs(wrapToPi(x)),30, 'Normalization', 'pdf')
k=pol(1)/Deff/2;
xfit=[0:pi/100:pi];
norm=1/2*sqrt(pi/k)*erf(sqrt(k)*pi);
hfit=1/norm*exp(-k*xfit.^2);
hold on;
plot(xfit,hfit)

%xlim([0,pi/4])
%%
%dxcropped = dx(:,1:end-1);
dxcropped = dx;
dxcropped(TS==0)=NaN;
xwrap2pi=wrapToPi(x);
Lat.a.auto_co_reinforcement(sign(dxcropped), Csimu, 20, 15)

%% check inter-bout distance

[binvals,elts_per_bin, v2bin, ~, ~ ] = ...
    BinsWithEqualNbofElements(dx, lambda, 20, 22 );
meanlambda = mean(v2bin,2);
varlambda = var(v2bin,1,2);
figure
plot(binvals,meanlambda);
hold on
plot(binvals, varlambda);

figure
histogram(lambda(TS==1), 'Normalization', 'pdf')
hold on
histogram(lambda(TS==0), 'Normalization', 'pdf')