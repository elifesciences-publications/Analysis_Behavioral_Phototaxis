%% Simulation sterovisual phototaxis

Nexp=10000; % number of independent experiments
Nsteps=50;  % number of time steps

w1=wturn;   % std of turning bouts amplitude distribution 
w2=wfor;    % std of forward bouts amplitude distribution 
p=pturn;    % probability of turning bout
pflip=0.21;  % neutral probability of flipping orientation
%boutn = 1:Nsteps;
%bias_slope = 20*(0.22-0.1*(log(-boutn/4))); % ratio of bias to contrast
%bias_slope(bias_slope<0) = 0;
bias_slope = 1.5;

%%
% initiate values 
lrst=(rand(Nexp,1)>0.5)*2-1; % initial left/right state
dx=[];  % change in orientation
x=-(rand(Nexp, 1)-0.5)*2*pi; % uniform initial orientation
qq=[]; % flipping probability
lr=[]; % left right state
% loop on tsteps
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
end

%% get corresponding contrast
percPm =1;
lum_lin = luminosity_linear(percPm);
L = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),pi-abs(pi-mod(x-pi/2,2*pi))); 
R = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),abs(pi-mod(x-pi/2,2*pi))); % experiment shifted by pi/2
Csimu = L-R;
Csimu = Csimu/max(abs(Csimu(:)));

%% polar histogram
if ~exist('XLat','var')
    path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
    name = 'lateralized_exps.mat';
    load([path name], 'El')
    XLat = El.AngleSource;
    Xfilt = El.AngleSourceFiltered;
    clear El
end

Nbins = 32;
angle = ((1:Nbins) - 0.5) * 2*pi/Nbins;
xpdf_data = hist(mod(XLat(:),2*pi),Nbins)/sum(hist(mod(XLat(:),2*pi),Nbins));
XLsub = XLat([250:end],2:22);
xpdf_datasub = hist(mod(XLsub(:),2*pi),Nbins)/sum(hist(mod(XLsub(:),2*pi),Nbins));
xpdf_simu = hist(mod(x(:)-pi/2,2*pi),Nbins)/sum(hist(mod(x(:)-pi/2,2*pi),Nbins));

%***
figure
polarplot(angle, xpdf_datasub, 'Linewidth', 1.5, 'DisplayName', 'subdata')
hold on
polarplot(angle,xpdf_data, 'Linewidth', 1.5, 'DisplayName', 'data')
polarplot(angle,xpdf_simu, 'Linewidth', 1.5, 'DisplayName', 'simu')
legend

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
xwrap2pi=wrapToPi(x(:,1:end-1));
%[Xbin,Drift, ~, ~ ] = histprofile(stereo_bias(xwrap2pi(:)), dx(:),20 );
[Xbin,elts_per_bin, v2bin, ~, ~ ] = ...
    BinsWithEqualNbofElements(Csimu(:,1:end-1), dx, 10, 12 );
Drift = mean(v2bin,2);
pol = polyfit(Xbin,Drift,1)
Drift_fit= polyval(pol,Xbin);
figure
plot(Xbin,Drift);
hold on;
plot(Xbin, Drift_fit)
%DD=[DD, p(1)]
%end
pflip*bias_slope*sqrt(2/pi)*w1*p/pol(1)
%ylim([-0.1,0.1])
hold on;
plot([-1 1], [0 0])
plot([0 0] , [-.1 .1])

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
xwrap2pi=wrapToPi(x);
Lat.auto_co_reinforcement(dx, Csimu, pturn, wturn, 20, 15)