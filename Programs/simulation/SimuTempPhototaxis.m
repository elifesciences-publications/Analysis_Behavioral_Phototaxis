
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

[dx_histogram, x_histogram] = hist(alldx, 100);
f = fit(x_histogram.',dx_histogram.','gauss2'); 

%***
fig = figure;
plot(x_histogram, dx_histogram, 'Linewidth', 2)
hold on
forward = @(x) f.a1*exp(-((x-f.b1)/f.c1).^2);
side = @(x) f.a2*exp(-((x-f.b2)/f.c2).^2);
plot(x_histogram,  forward(x_histogram)+side(x_histogram), x_histogram, side(x_histogram), 'Linewidth', 1.5)

%% check proba functions
dll = -2:0.01 : 1;
for i = 1 : length(dll)
    [pt(i), wt(i)] = ProbaTurnFLum(dll(i));
end

%***
figure;
plot(dll, pt, dll, wt);
legend('proba turn', 'w turn')

%% lum
lum_exp6 = luminosity_exponentielle(0.6);
lum_exp3 = luminosity_exponentielle(0.3);
lum_sin6 = luminosity_sinus(0.6);
lum = lum_exp6;

% intitialization
%--------------------------------------------------------------------------
psw_side = 0.188; % probabibility of switching left vs right states
p_turn = 0.4; % probability of triggering a turn swim (otherwise go straight)

wturn=0.69; % 0.74
wstraight=0.1; % 0.12
%--------------------------------------------------------------------------

%rng('shuffle')
Nexp = 10000; % number of experiments
Ntimes = 50; % number of time steps per experiment
theta_complete = NaN(Nexp, Ntimes); % complete dataset
luminosity_complete = NaN(Nexp, Ntimes);

theta_ini = rand(1,Nexp)*2*pi; % initial angles uniformly distributed
luminosity_ini = interp1(lum(:,1), lum(:,2), abs(wrapTo180(rad2deg(theta_ini))));
lrst = (rand > 0.5)*2-1; % intial turn state random 50/50

tic
for N = 1:Nexp
    theta = NaN(1,Ntimes);
    luminosity = NaN(1,Ntimes);
    theta(1) = theta_ini(N);
    luminosity(1) = luminosity_ini(N);
    for t = 1:Ntimes-1
        th = theta(t);
        lrst = lrst*((rand > psw_side)*2-1); % independent from whole field illum
        if t > 2
            dll = (luminosity(t) - luminosity(t-1))/((luminosity(t) + luminosity(t-1))/2) ; 
            [p_turn, wturn] = ProbaTurnFLum(dll);
        end
        turn = rand < p_turn;
        if turn
            dth = lrst * abs(wturn*randn);
        else
            dth = (wstraight*randn);
        end
        
        theta(t+1) = th + dth;
        
        luminosity(t+1) = interp1(lum(:,1), lum(:,2), abs(wrapTo180(rad2deg(theta(t+1)))));
    end
    theta_complete(N,:) = theta;
    luminosity_complete(N,:) = luminosity;
end

alpha = theta_complete(:);
[mp, rho_p, mu_p] = circ_moment(alpha);


Nbins = 18;
angle = ( (1:Nbins) -0.5) * 2*pi/Nbins;
PDFangle = hist(mod(theta_complete(:),2*pi),Nbins)/sum(hist(mod(theta_complete(:),2*pi),Nbins));
PDFangleExp =  hist(mod(e6_x(:),2*pi),Nbins)/sum(hist(mod(e6_x(:),2*pi),Nbins));
%PDFangleExp =  hist(mod(e3_x(:),2*pi),Nbins)/sum(hist(mod(e3_x(:),2*pi),Nbins));
%PDFangleExp =  hist(mod(s6_x(:),2*pi),Nbins)/sum(hist(mod(s6_x(:),2*pi),Nbins));

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
med_bouts_per_seq = median(size(e6_x,2)-sum(isnan(e6_x),2))+10;
xData = e6_x(:,1:med_bouts_per_seq);
xSimu = theta_complete;
f = e6_f;

colourID = 2;
[fig] = polar_distribution_simu_exp(xData+pi, xSimu+pi, Nbins, f, colourID)


toc
%% time evolution
b = 1 : 10 : Ntimes;
for i = 2 : Ntimes
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
