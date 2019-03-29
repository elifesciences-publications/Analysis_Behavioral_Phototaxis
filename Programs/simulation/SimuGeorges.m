
%% ORIENTATION
%%
rng('shuffle')
Nexp = 100; % number of experiments
Ntimes = 1000; % number of time steps per experiment
theta_complete = NaN(Nexp, Ntimes); % complete dataset
theta_ini = rand(1,Nexp)*360; % initial angles for all experiments
%theta_ini=-90*ones(1,Nexp);
bias = @(x) 2*((mod( (x), 360) > 180-1 ) >0 ) - 1; % sign of the bias (towards 0)
DP = 0;%-0.1; % difference between p+ and p-
%beat_amp = 33;

% sin profile
sigba = 25; 
beat_ampSin60 = @(x) 25 + 33*exp( -(x/sigba).^2) ;

% exp profile
sigba = 75; 
beat_ampExp60 = @(x) 20 + 20*exp(- (x/sigba).^2 );

dth = NaN(Nexp, Ntimes);
tic
randthresh = 0.5;%[0.4: 0.005 : 0.6];
for i = 1 : length(randthresh)
for N = 1:Nexp
    theta = theta_ini(N);
    for t = 1:Ntimes-1
        th = theta(t);
        if DP ~= 0
            step = 2*double( rand(1) > (DP/2*bias(th)+0.5) )-1;
        else
            step = 2*double( rand >  randthresh(i)) - 1;
        end
        noise = 1+0.3*randn ; % I add noise to avoid artifact
        dth(N,t) = step*beat_ampExp60(wrapTo180(th))*noise;
        theta(t+1) = th + dth(N,t);  
    end
    theta_complete(N,:) = theta;
end
alpha = theta_complete(:)/180*pi;
[mp, rho_p, mu_p(i)] = circ_moment(alpha);
ddth(i) = sum(dth(:)>0) - sum(dth(:)<0);
end

Nbins = 36;
angle = ( (1:Nbins) -0.5) * 360/Nbins;
PDFangle = hist(mod(theta_complete(:),360),Nbins)/sum(hist(mod(theta_complete(:),360),Nbins));
%PDFangle=hist(mod(theta_comp(:),360),angle)

%fig = figure;
subplot(2,1,1)
plot(angle,PDFangle);
subplot(2,1,2)
polarplot([angle/180*pi,angle(1)/180*pi],[PDFangle,PDFangle(1)]);
set(gca,'ThetaZeroLocation','left',...
        'ThetaDir','clockwise')
    title('sin')

%orient_mean=cos(angle/180*pi)*PDFangle'
%polarplot(hist(mod(theta_comp(:),360),19))
%plot(cumsum(hist(mod(theta_comp(:),360),19)/sum(hist(mod(theta_comp(:),360),19))))

% calculate mean angle and moment of the distribution.
display(['mean angle in degree  :  ',num2str(mu_p*180/pi)]);
display(['first moment  :  ',num2str(rho_p)]);
[p,z] = circ_rtest(theta_complete(:)/180*pi);
disp(['p = ' num2str(p) ' ; z = ' num2str(z)])
 
toc
%% LUMINANCE

rng('shuffle')
Nexp = 1000; % number of experiments
Ntimes = 50; % number of time steps per experiment
theta_complete = NaN(Nexp, Ntimes); % complete dataset
theta_ini = rand(1,Nexp)*360; % initial angles for all experiments

% bout amplitude from luminance
lum_exp6 = luminosity_exponentielle(0.6);
                
sigba = 0.02;
ba = @(lu) 0.1 + 1.4*exp(- (lu/sigba).^2 );

dth = NaN(Nexp, Ntimes);
tic
for N = 1:Nexp
    theta = theta_ini(N);
    for t = 1:Ntimes-1
        th = abs(wrapTo180(theta(t)));
        luminosity = interp1(lum_exp6(:,1), lum_exp6(:,2), th);
        dth(N,t) = 2*double( rand >  randthresh(i)) - 1;
        noise = 1+0.3*randn ;
        theta(t+1) = th + ba(luminosity)*noise*dth(N,t); 
    end
    theta_complete(N,:) = theta;
end
alpha = theta_complete(:)/180*pi;
[mp, rho_p, mu_p(i)] = circ_moment(alpha);
ddth(i) = sum(dth(:)>0) - sum(dth(:)<0);


Nbins = 36;
angle = ( (1:Nbins) -0.5) * 360/Nbins;
PDFangle = hist(mod(theta_complete(:),360),Nbins)/sum(hist(mod(theta_complete(:),360),Nbins));
%PDFangle=hist(mod(theta_comp(:),360),angle)

fig = figure;
subplot(2,1,1)
plot(angle,PDFangle);
subplot(2,1,2)
polarplot([angle/180*pi,angle(1)/180*pi],[PDFangle,PDFangle(1)]);
set(gca,'ThetaZeroLocation','left',...
        'ThetaDir','clockwise')
    title('sin')

%% calculate drift
drift_comp=zeros(Nexp, Ntimes-1);
%times=[1:Ntimes-1];
for N=1:Nexp
    theta=theta_complete(N,:);
    dth=diff(theta);
    bias_light=bias(theta(1:end-1));
    drift_comp(N,:)=cumsum(dth.*bias_light);
   
   %plot(drift)
   %drift/times
   %DP*beat_amp
    
end
 
figure
plot(mean(drift_comp,1));
 
%% analyse data
%addpath('/Users/georges/Documents/matlab/CircStat2012a')
Nbins=18;
angle=([1:Nbins]-0.5)*360/Nbins;
hold off;
PDFangle=hist(ang_lab_sin_ok(:),angle)/sum(hist(ang_lab_sin_ok,Nbins));
polarplot([angle/180*pi,angle(1)/180*pi],[PDFangle,PDFangle(1)]);
 
hold on;
PDFangle=hist(ang_lum_sin(:),Nbins)/sum(hist(ang_lum_sin,Nbins));
polarplot([angle/180*pi,angle(1)/180*pi],[PDFangle,PDFangle(1)]);

