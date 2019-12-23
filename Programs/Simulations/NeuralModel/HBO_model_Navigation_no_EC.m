%clear all

%% --- Load Data ---
% load('/Users/sebastienwolf/Documents/Project/Neurofish/Programs/Sebastien/HBO 3.0/BRAIN region model/inter_bout.mat')
% load('/Users/sebastienwolf/Documents/Project/Neurofish/Programs/Sebastien/HBO 3.0/BRAIN region model/dist.mat')
inter_bout_exp = IBIi(:);
inter_bout_exp(isnan(inter_bout_exp))=[];

dist = Disti;
dist(isnan(dist)==1)=[];

io_in_list = [0:1000:25000];%100000];%1000;
WE_list = [0.7:0.01:1];
C_list = [-1:0.1:1];
Itot_list = [500 750 1000];

%% --- Define parameters ---
%Parameters you can define
Param.T = 1000000 %duration  of the simulation in bins
Param.sigma = 500 %noise of the model
Param.io_in = 0
Param.dt = 0.1 %frame rate in s
Param.sig_scout = 0.1
Param.sig_turn = 0.6
Param.pturn = 0.41
Param.WE = 0.962; %0.97
Param.WI = 7
Param.tau = 0.1
Param.i0 = 20
Param.Tporte = 0.1
Param.Itot = 10000; %1100 %Itot_list(p_ilist)
Param.teta0 = pi; %origin
Param.x_fish0 = 0; %origin
Param.y_fish0= 0; %origin


% load parameters and initialization
io_in = Param.io_in;
T = Param.T;
sigma = Param.sigma;
dt = Param.dt;
sig_scout = Param.sig_scout;
sig_turn = Param.sig_turn;
WE = Param.WE;
WI = Param.WI;
tau = Param.tau;
io = Param.i0;
Tporte = Param.Tporte;
increm_C = 0;
teta = Param.teta0;
teta_cum = teta;
C = 2/pi*(abs(pi/2-teta)-pi/2);
Itot = Param.Itot;
Light_cur_R = Itot*(1+C)/2;
Light_cur_L = Itot*(1-C)/2;


%% --- Simulation ---
%sampling in the distribution of interbouts
spikeL = zeros(1,T);
spike_forward = zeros(1,T);
spikeL(round(1/dt)) = 1;
spikeL_all_bout = spikeL;
pos_bout = round(1/dt);
s_bout = 1;

% find position of bouts pos_bout in time
while   pos_bout<T-max(inter_bout_exp)/dt
    pos_bout = pos_bout+round(inter_bout_exp(randi(length(inter_bout_exp)))/dt);
    if pos_bout<T
        s_bout = s_bout+1;
        type_spikeL_all_bout(s_bout) = double(rand()<Param.pturn);
        spikeL_all_bout(pos_bout) = 1;
        
        time_bout(s_bout) = pos_bout;
        
        if type_spikeL_all_bout(s_bout)==1
            spikeL(pos_bout) = 1;
        end
        
        if type_spikeL_all_bout(s_bout)==0
            spike_forward(pos_bout) = 1;
        end
    end
end
size(type_spikeL_all_bout)

stimL = find(spikeL>0);
stimL_all_bout = find(spikeL_all_bout>0);
length(stimL_all_bout)


for t= 1:round(Tporte/dt)
    t;
    spikeL(find(spikeL==1)+t)=1;
end

spikeL(T+1:end)=[];
spikeR=spikeL;

%initialisation
isl=zeros(1,T);
isr=zeros(1,T);
rr=zeros(1,T);
rl=zeros(1,T);

s_bout=0;
for i=1:T-1 
    if rl(i)-rr(i)>0
        spikeR(i)=0;
    else
        spikeL(i)=0;
    end
    
    %rl
    xi=sigma*randn(1);
    arg=WE*rl(i)-WI*rr(i)+io+isl(i)+io_in*spikeL(i)+Light_cur_L(i);
    rl(i+1)=rl(i)+dt*(1/tau)*(-rl(i)+heaviside(arg)*arg+xi);
    
    %rr
    xi=sigma*randn(1);
    arg=WE*rr(i)-WI*rl(i)+io+isr(i)+io_in*spikeR(i)+Light_cur_R(i);
    rr(i+1)=rr(i)+dt*(1/tau)*(-rr(i)+heaviside(arg)*arg+xi);
    
    % left side -
    if rl(i)-rr(i)>0 
        
        if spikeL(i)>0
            dteta = sig_turn*abs(randn(1));
            teta(i+1) = wrapToPi(teta(i) - dteta - pi/2)+pi/2 ;
            teta_cum(i+1) = teta_cum(i) - dteta ;
            s_bout = s_bout+1;
            teta_bout(s_bout) = teta(i+1);
            teta_cum_bout(s_bout) = teta_cum(i+1);
            bout_type(s_bout) = -1;
            bout_time(s_bout) = i*dt;
        end
        if spike_forward(i)>0
            dteta = sig_scout*abs(randn(1));
            teta(i+1) = wrapToPi(teta(i) - dteta - pi/2)+pi/2;
            teta_cum(i+1) = teta_cum(i) - dteta ;
            s_bout = s_bout+1;
            teta_bout(s_bout) = teta(i+1);
            teta_cum_bout(s_bout) = teta_cum(i+1);
            bout_type(s_bout) = 0;
            bout_time(s_bout) = i*dt;
        end
        if spikeL(i)==0 && spike_forward(i)==0
            teta(i+1) = teta(i);
            teta_cum(i+1) = teta_cum(i);
        end
    end
    
    % right side +
    if rr(i)-rl(i)>0
        
        if spikeR(i)>0
            dteta = sig_turn*abs(randn(1));
            teta(i+1) = wrapToPi(teta(i) + dteta - pi/2)+pi/2 ;
            teta_cum(i+1) = teta_cum(i) + dteta ;
            s_bout = s_bout+1;
            teta_bout(s_bout) = teta(i+1);
            teta_cum_bout(s_bout) = teta_cum(i+1);
            bout_type(s_bout) = 1;
            bout_time(s_bout) = i*dt;
        end
        if spike_forward(i)>0
            dteta = sig_scout*abs(randn(1));
            teta(i+1) = wrapToPi(teta(i) + dteta - pi/2)+pi/2;
            teta_cum(i+1) = teta_cum(i) + dteta ;
            s_bout = s_bout+1;
            teta_bout(s_bout) = teta(i+1);
            teta_cum_bout(s_bout) = teta_cum(i+1);
            bout_type(s_bout) = 0;
            bout_time(s_bout) = i*dt;
        end
        if spikeR(i)==0 && spike_forward(i)==0
            teta(i+1) = teta(i);
            teta_cum(i+1) = teta_cum(i);
        end
    end
    
    
    if rl(i)==rr(i)
        teta(i+1) = teta(i);
        teta_cum(i+1) = teta_cum(i);
    end
    
    C(i+1)=2/pi*(abs(pi/2-teta(i+1))-pi/2);
    
    Light_cur_R(i+1) = Itot*(1+C(i+1))/2;
    Light_cur_L(i+1) = Itot*(1-C(i+1))/2;
    
    if C(i+1)>=1; C(i+1)=1;end
    if C(i+1)<=-1; C(i+1)=-1;end
    
end

x_fish = zeros(1,length(type_spikeL_all_bout));
y_fish = zeros(1,length(type_spikeL_all_bout));
x_fish(1) = Param.x_fish0;
y_fish(1) = Param.y_fish0;

for s_bout = 1:length(type_spikeL_all_bout)-1
    dist_bout(s_bout) = dist(randi(length(dist)));
    x_fish(s_bout+1) = x_fish(s_bout) + dist_bout(s_bout)*cos(teta_bout(s_bout));
    y_fish(s_bout+1) = y_fish(s_bout) + dist_bout(s_bout)*sin(teta_bout(s_bout));
end

%% plots of interest

figure
times=0:dt:dt*(T-dt);
hold on
plot(times,10000*spikeL,'r:','LineWidth',2)
plot(times,10000*spikeR,'b:','LineWidth',2)
plot(times,rl,'LineWidth',2,'color','r')
plot(times,rr,'LineWidth',2,'color','b')

figure;
plot(x_fish*0.1,y_fish*0.1,'k ')
hold on
plot(x_fish(1:1000)*0.1,y_fish(1:1000)*0.1,'r ')

% show raw signal
%***
figure
hold on
plot(times,rl,'LineWidth',2,'color','r')
plot(times,rr,'LineWidth',2,'color','b')
title('raw')
xlim([0 50])

% transform into telegraph signal
rlbin = rl;
rrbin = rr;
rlbin(rl>rr) = 1;
rlbin(rl<rr) = 0;
rrbin(rr>rl) = 1;
rrbin(rr<rl) = 0;
telsig = rlbin - rrbin;

%***
figure
plot(times, telsig)
title('telgraph signal')
xlim([0 50])
ylim([-1.5 1.5])

% find residence time
ds = diff(telsig);
stateonset = find(ds); % loc of flips
restime = [stateonset(1) diff(stateonset)]; % temps de résidence
restime = restime*dt; % en secondes 
restimex = [0:0.1:80];
restimepdf = computePDF(restimex, restime);

%***
figure
plot(restimex, restimepdf)
title('residence time of telegraph signal')
ylabel('pdf')

%***
turn_thresh = 0.22;
db = 10;

diffteta = diff(teta_cum_bout);
dtb = diff(bout_time);
turns = spikeL + spikeR;
[fig1, fig2, binvals1, mv1, stdv1] = Spont.autocorrelationVSinterboutinterval_nonconsecutive...
    (diffteta, dtb, turn_thresh, db, 'simu seb ' );

% --- bias with contrast ---
dX = diff(teta_cum_bout);
contrast = C(isfinite(bout_type));
Vart1 = contrast(2:end);
Vart2 = dX;
b=10;

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, b, b+3);
mV2 = nanmean(v2bin,2);
sV2sq = nanstd(v2bin,1,2);

errorbar(binvals, mV2, sV2sq/sqrt(elts_per_bin-1))

% --- pflip with contrast ---
dX = diff(teta_cum_bout);
contrast = C(isfinite(bout_type));

binsdXdX = 9;
binspflip = 7;
Lat.auto_co_reinforcement(dX, dX, contrast, binsdXdX, binspflip, 0, 'absolute')
