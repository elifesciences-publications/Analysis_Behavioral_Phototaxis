
% load data
path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
name = 'lateralized_exps.mat';
load([path name], 'El')
XLat = El.AngleSource;

%%
psw_turn = 0.2; % probabibility of switching left vs right states
p_turn = 0.55; % probability of triggering a turn swim (otherwise go straight)

wturn= 0.74; % 1.1
wstraight= 0.11; % 0.13

mturn = @(x) -pi/12*x; % pi/12
%%
% intitialization
% -------------------------------------------------------------------------
rng('shuffle')
Nexp = 10000; % number of experiments
Ntimes = 30; % number of time steps per experiment
theta_complete = NaN(Nexp, Ntimes); % complete dataset
lL_R = NaN(Nexp, Ntimes);

initial_bias = 0;
if logical(initial_bias)
    theta_ini = randsample(XLat(:,1), Nexp, true); % bias the initial distribution
else
    theta_ini = rand(1,Nexp)*2*pi; % initial angles uniformly distributed
end
lum_lr_ini = interp1([-pi 0 pi], [1 -1 1], wrapToPi(theta_ini));
%lrst = (rand > 0.5)*2-1; % intial turn state random 50/50
lrst = 1;

tic
for N = 1:Nexp
    theta = NaN(1,Ntimes);
    lum_lr = NaN(1,Ntimes);
    theta(1) = theta_ini(N);
    lum_lr(1) = lum_lr_ini(N);
    for t = 1:Ntimes-1
        th = theta(t);
        %lrst = lrst*((rand > psw_turn)*2-1); % independant from whole field illum
        turn = rand <= p_turn;
        if turn
            dth = lrst * (wturn*randn) - mturn(lum_lr(t));
        else
            dth = lrst * (wstraight*randn);
        end
        
        theta(t+1) = th + dth;
        
        lum_lr(t+1) = interp1([0 pi], [-1 1], abs(wrapToPi(theta(t+1))));
    end
    theta_complete(N,:) = theta;
    lL_R(N,:) = lum_lr;
end

THETA = theta_complete(:);
[mp, rho_p, mu_p] = circ_moment(THETA);

% ***
figure
Nbins = 32;
angle = ( (1:Nbins) -0.5) * 2*pi/Nbins;
PDFangle = hist(mod(theta_complete(:),2*pi),Nbins)/sum(hist(mod(theta_complete(:),2*pi),Nbins));
Xexp = XLat(:, 1:20);
PDFangleLat =  hist(mod(Xexp(:),2*pi),Nbins)/sum(hist(mod(Xexp(:),2*pi),Nbins));

% ***
fig = figure;

polarplot([angle,angle(1)]-pi/2,[PDFangle,PDFangle(1)], 'k', 'Linewidth', 2);
hold on
polarplot([angle,angle(1)]-pi/2,[PDFangleLat,PDFangleLat(1)], 'Linewidth', 1);
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','clockwise')
legend('simu', 'exp')
title('Lateralized : 10000 sequences, 30 bouts')
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
[p,z] = circ_rtest(theta_complete(:)/180*pi);
disp(['p = ' num2str(p) ' ; z = ' num2str(z)])
 
toc

meanX = circ_mean(theta_complete);
R = circ_r(theta_complete);
Rproj_perfish = R.*cos(meanX+pi/2);

figure; 
plot(R)
hold on
plot(Rproj_perfish)

disp('R proj')
circ_r(theta_complete(:))*cos(circ_mean(theta_complete(:)+pi/2)) 

%% bout bias
thetaForDiff = wrapToPi(theta_complete(:,1:end-1));
dtheta = diff(theta_complete, 1, 2 );
lL_Rfd = lL_R(:,1:end-1);

% on orientation
wo = pi/12;
bins_orient = wo/2 : wo : pi ;
bins_orient = [sort(-bins_orient), bins_orient];
meandX_orient = NaN(length(bins_orient), 1);
semdX_orient = NaN(length(bins_orient), 1);
vardX_orient = NaN(length(bins_orient), 1);
for i = 2 : length(bins_orient)
    dXselected = dtheta(thetaForDiff < bins_orient(i) & thetaForDiff > bins_orient(i-1));
    n = length(dXselected);
    semdX_orient(i-1) = nanstd(dXselected)/sqrt(n);
    meandX_orient(i-1) = nanmean(dXselected);
    vardX_orient(i-1) = nanvar(dXselected);
end

% on luminosity
wi = 0.05;
bins_ill = -1 : wi : 1 ;
meandX_ill = NaN(length(bins_ill), 1);
semdX_ill = NaN(length(bins_ill), 1);
vardX_ill = NaN(length(bins_ill), 1);
for i = 2 : length(bins_ill)
    dXselected = dtheta(lL_Rfd < bins_ill(i) & lL_Rfd > bins_ill(i-1));
    n = length(dXselected);
    semdX_ill(i-1) = nanstd(dXselected)/sqrt(n);
    meandX_ill(i-1) = nanmean(dXselected);
    vardX_ill(i-1) = nanvar(dXselected);
end

%***
f = figure;
subplot(1,2,1)
yyaxis left
%scatter(thetaForDiff(:), dtheta(:))
hold on
plot((bins_orient + wo/2), (meandX_orient), 'Linewidth', 2)
errorbar((bins_orient + wo/2), (meandX_orient), (semdX_orient), '.', 'Color', [0.2 0.2 0.2])
ylabel('dX')

yyaxis right
pvo = plot((bins_orient + wo/2), (vardX_orient));
pvo.Color(4) = 0.3;
ylabel('dX^2')
grid on

subplot(1,2,2)
%scatter(lL_Rfd(:), dtheta(:));
hold on
plot((bins_ill + wi/2), (meandX_ill), 'Linewidth', 2)
errorbar((bins_ill + wi/2), (meandX_ill), (semdX_ill), '.', 'Color', [0.2 0.2 0.2])
ylabel('dX')

yyaxis right
pvi = plot((bins_ill + wi/2), (vardX_ill), 'Linewidth', 1);
pvi.Color(4) = 0.3;

subplot(1,2,1)
xlabel('X (rad)')
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi' '-\pi/2' '0' '\pi/2' '\pi'})
subplot(1,2,2)
xlabel('Lum(L)-Lum(R)')