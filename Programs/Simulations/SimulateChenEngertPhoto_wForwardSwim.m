%% Simulation sterovisual phototaxis

wturn=0.74; % 0.74
wfor=0.11; % 0.12
pflip=0.18;  % neutral probability of flipping orientation
pturn = 0.4;
bias_slope = 1.5;

%% check proba functions
dll = -2:0.01 : 1;
[pt, wt] = ProbaTurnFLum(dll);

%***
figure;
plot(dll, pt, dll, wt);
legend('proba turn', 'w turn')
%% generate two populations of inter-bout distances
sigmafwd = 0.6;
sigmaturn = 0.5;
muturn=1;
mufwd = 0.7;
possible_lambda_values = 0:0.01:30;

%%
Nexp=50; % number of independent experiments
Nsteps=10000;  % number of time steps
virtual_roi_size = 50 ;

% initiate values 
dTheta = NaN(Nexp, Nsteps);  % change in orientation
Theta = NaN(Nexp, Nsteps);
Theta(:,1) = -(rand(Nexp, 1)-0.5)*2*pi; % uniform initial orientation
Q = NaN(Nexp, Nsteps); % flipping probability
LRstate = NaN(Nexp, Nsteps); % left right state
Lambda = NaN(Nexp, Nsteps);
TS = zeros(Nexp, Nsteps);
Lum = ones(Nexp, Nsteps);
Dll = NaN(Nexp, Nsteps);
a = zeros(Nexp, Nsteps);
b = zeros(Nexp, Nsteps);

lrst = (rand(Nexp,1)>0.5)*2-1;
pturn = 0.4;
% loop on tsteps
for t = 1 : Nsteps
    % --- pflip : stays constant in tempo ---
    q = pflip; %*(1-0.5*bias_slope*lrst.*stereo_bias_mat(thetalast));
    q(q<0) = 0;
    q(q>1) = 1;
    lrst = lrst.*((rand(Nexp,1)>q)*2-1);
    % --- variables determining turn angle : pturn, wturn ---
    if t > 2
        dll = (Lum(:,t) - Lum(:,t-1))./((Lum(:,t) + Lum(:,t-1))/2) ;
    else
        dll = zeros(Nexp,1);
    end
    [pturn, wturn] = ProbaTurnFLum(dll);
    Dll(:,t) = dll;
    turn = rand(Nexp,1) < pturn;
    indturn = find(turn==1);
    indfor = find(turn==0);
    % --- turn angle ---
    dtheta = zeros(Nexp,1);
    dtheta(indturn) = lrst(indturn).*abs(wturn(indturn).*randn(length(indturn),1));
    dtheta(indfor) = wfor*randn(length(indfor),1);
    % --- inter-bout distance ---
    l = NaN(Nexp,1);
    l(indturn) = pick_from_lognorm_distrib(possible_lambda_values, sigmaturn, muturn, length(indturn))';
    l(indfor) = pick_from_lognorm_distrib(possible_lambda_values, sigmafwd, mufwd, length(indfor))';
    % --- cartesian coordinates ---
    if t > 1
        [da, db] = pol2cart(Theta(:,t-1), l);
        a(:,t) = a(:,t-1) + da;
        b(:,t) = b(:,t-1) + db;
        Theta(:,t) = Theta(:,t-1) + dtheta;
    end
    out_of_roi = (sqrt(a(:,t).^2 + b(:,t).^2) > virtual_roi_size/2);
    Lum(out_of_roi, t) = 0;
    
    % store variables
    dTheta(:,t) = dtheta;
    LRstate(:,t) = lrst;
    Lambda(:,t) = l;
    TS(:,t) = turn;
end


%% get corresponding contrast
percPm =1;
lum_lin = luminosity_linear(percPm);
L = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),pi-abs(pi-mod(Theta-pi/2,2*pi))); 
R = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),abs(pi-mod(Theta-pi/2,2*pi))); % experiment shifted by pi/2
Csimu = L-R;
Csimu = Csimu/max(abs(Csimu(:)));

%%


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
xpdf_simu = hist(mod(Theta(:)-pi/2,2*pi),Nbins)/sum(hist(mod(Theta(:)-pi/2,2*pi),Nbins));

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
xpart=abs(wrapToPi(Theta(:,tinit:tend)));
histogram(xpart(:),15, 'Normalization', 'pdf');

%%

Nexp_real=300; % number of trajectories in actual experiments;
for pturn=1:1000;
    xpart=Theta(fix(Nexp*abs(rand(Nexp_real,1)))+1,tinit:tend);
    hval(pturn,:)=histcounts(abs(wrapToPi(xpart(:))),binedges, 'Normalization', 'pdf');
    %hval(p,:)=hpart.Values;
end
hstd=std(hval,1);
figure;
boundedline([1:15]/15*pi,hmean_val, hstd);
%% extract the drift
xwrap2pi=wrapToPi(Theta(:,1:end-1));
%[Xbin,Drift, ~, ~ ] = histprofile(stereo_bias(xwrap2pi(:)), dx(:),20 );
[binvals,elts_per_bin, v2bin, ~, ~ ] = ...
    BinsWithEqualNbofElements(Csimu(:,1:end-1), dTheta, 10, 12 );
Drift = mean(v2bin,2);
pol = polyfit(binvals,Drift,1)
Drift_fit= polyval(pol,binvals);
figure
plot(binvals,Drift);
hold on;
plot(binvals, Drift_fit)
%DD=[DD, p(1)]
%end
pflip*bias_slope*sqrt(2/pi)*wturn*pturn/pol(1)
%ylim([-0.1,0.1])
hold on;
plot([-1 1], [0 0])
plot([0 0] , [-.1 .1])

%% extract the drift
dtheta = dTheta;
dtheta(TS==0) = NaN;
dtheta(dtheta<0) = -1; dtheta(dtheta>0) = 1;
ac1 = dtheta(:,1:end-1).*dtheta(:,2:end)
[binvals,elts_per_bin, v2bin, ~, ~ ] = ...
    BinsWithEqualNbofElements(Csimu(:,1:end-2), ac1, 10, 12 );
mac1 = mean(v2bin,2);
figure;
plot(binvals, mac1)

%% get effective diffusion coefficient

dxbias=pol(1)*wrapToPi(Theta(:,1:end-1));
plot(wrapToPi(Theta(:,1:end-1)), dxbias);
xdiff=cumsum(dTheta+dxbias);
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
histogram(abs(wrapToPi(Theta)),30, 'Normalization', 'pdf')
k=pol(1)/Deff/2;
xfit=[0:pi/100:pi];
norm=1/2*sqrt(pi/k)*erf(sqrt(k)*pi);
hfit=1/norm*exp(-k*xfit.^2);
hold on;
plot(xfit,hfit)

%xlim([0,pi/4])
%%
%dxcropped = dx(:,1:end-1);
dxcropped = dTheta;
dxcropped(TS==0)=NaN;
xwrap2pi=wrapToPi(Theta);
Lat.a.auto_co_reinforcement(sign(dxcropped), Csimu, 20, 15)

%% check inter-bout distance

[binvals,elts_per_bin, v2bin, ~, ~ ] = ...
    BinsWithEqualNbofElements(dTheta, Lambda, 20, 22 );
meanlambda = mean(v2bin,2);
varlambda = var(v2bin,1,2);
figure
plot(binvals,meanlambda);
hold on
plot(binvals, varlambda);

figure
histogram(Lambda(TS==1), 'Normalization', 'pdf')
hold on
histogram(Lambda(TS==0), 'Normalization', 'pdf')