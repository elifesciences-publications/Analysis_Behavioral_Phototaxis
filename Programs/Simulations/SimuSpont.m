% simulation unbiased navigation

%%
psw_turn = 0.1; % probabibility of switching left vs right states
p_turn = 0.5; % probability of triggering a turn swim (otherwise go straight)

wturn=0.74; % 0.74
wstraight=0.12; % 0.12

%%
% intitialization
% -------------------------------------------------------------------------
rng('shuffle')
Nexp = 10000; % number of experiments
Ntimes = 100; % number of time steps per experiment
theta_complete = NaN(Nexp, Ntimes); % complete dataset
luminosity_complete = NaN(Nexp, Ntimes);

initial_bias = 1;
start = 5;
if initial_bias ==1
    theta_ini = randsample(e6_x(:,start), Nexp, true); % bias the initial distribution
else
    theta_ini = rand(1,Nexp)*2*pi; % initial angles uniformly distributed
end
lrst = (rand > 0.5)*2-1; % intial turn state random 50/50

tic
for N = 1:Nexp
    theta = NaN(1,Ntimes);
    theta(1) = theta_ini(N);
    for t = 1:Ntimes-1
        th = theta(t);
        lrst = lrst*((rand > psw_turn)*2-1); % independant from whole field illum
        turn = rand < p_turn;
        if turn
            dth = lrst * abs(wturn*randn);
        else
            dth = lrst * abs(wstraight*randn);
        end
        theta(t+1) = th + dth;
    end
    theta_complete(N,:) = theta;
end

alpha = theta_complete(:);
[mp, rho_p, mu_p] = circ_moment(alpha);

% ***
figure
Nbins = 64;
angle = ( (1:Nbins) -0.5) * 2*pi/Nbins;
PDFangle = hist(mod(theta_complete(:),2*pi),Nbins)/sum(hist(mod(theta_complete(:),2*pi),Nbins));
PDFangleExp =  hist(mod(e6_x(:),2*pi),Nbins)/sum(hist(mod(e6_x(:),2*pi),Nbins));

%fig = figure;
subplot(2,1,1)
plot(angle,PDFangle);
hold on
plot(angle,PDFangleExp);
legend('simu', 'exp60')

subplot(2,1,2)
polarplot([angle,angle(1)],[PDFangle,PDFangle(1)]);
set(gca,'ThetaZeroLocation','left',...
        'ThetaDir','clockwise')
hold on
polarplot([angle,angle(1)],[PDFangleExp,PDFangleExp(1)]);
set(gca,'ThetaZeroLocation','left',...
        'ThetaDir','clockwise')

% calculate mean angle and moment of the distribution.
display(['mean angle in degree  :  ',num2str(mu_p*180/pi)]);
display(['first moment  :  ',num2str(rho_p)]);
[p,z] = circ_rtest(theta_complete(:)/180*pi);
disp(['p = ' num2str(p) ' ; z = ' num2str(z)])
 
toc

%% time evolution
% b = 1 : 10 : Ntimes;
% for i = 2 : Ntimes
%     polarhistogram(theta_complete(:,i), 50)
%     drawnow
%     text(30, 15, num2str(i))
%     pause(0.2)
% end

Rsimu = circ_r(theta_complete, [], [], 1);
mean_theta = circ_mean(theta_complete, [], 1);
Rpsimu = Rsimu.*cos(mean_theta-pi);

figure
plot(Rsimu)
hold on
plot(Rpsimu)

%% comparison with experiment
% load experiment
loadPooledExperimentsCell

% time dyn script
TemporalDynamics

subplot(3,1,1)
plot(start : (Ntimes + start -1) , Rsimu, 'Linewidth', 1.5)

subplot(3,1,2)
plot(start : (Ntimes + start -1) , Rpsimu, 'Linewidth', 1.5)
