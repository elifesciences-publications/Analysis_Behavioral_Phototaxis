% Code from Sébastien

%clear all

%load('/Users/sebastienwolf/Documents/Project/Neurofish/Programs/Sebastien/HBO 3.0/BRAIN region model/inter_bout.mat')
%load('/Users/sebastienwolf/Documents/Project/Neurofish/Programs/Sebastien/HBO 3.0/BRAIN region model/dist.mat')
%dist(isnan(dist)==1)=[];
inter_bout_exp = IBIi(:);
inter_bout_exp(isnan(inter_bout_exp))=[];

dist = Disti;
dist(isnan(dist)==1)=[];

io_in_list = [0:1000:25000];%100000];%1000;
WE_list = [0.7:0.01:1];
C_list = [-1:0.1:1];
Itot_list = [500 750 1000];


%Parameters you can define
Param.T = 1200000        %duration  of the simulation in bins
Param.sigma = 500       %noise of the model
Param.io_in = 0
Param.dt = 0.1          %frame rate in s
Param.sig_scout = 0.1
Param.sig_turn = 0.6
Param.pturn = 0.41
Param.WE = 0.962
Param.WI = 7
Param.tau = 0.1
Param.i0 = 20
Param.Tporte = 0.1
Param.Contrast = 0.4;   %C_list(p_contrast)
Param.Itot = 1500;         %Itot_list(p_ilist) 1500



% load parameters
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

C = Param.Contrast;
mean_bias = [];

figure
hold on
%Loop on different contrasts
for c = [-0.5 0 0.5] 
    tic
    increm_C = 0;
    Itot = Param.Itot;
    Light_cur_R = Itot*(1+c)/2;
    Light_cur_L = Itot*(1-c)/2;
    
    %sampling in the distribution of interbouts
    spikeL = zeros(1,T);
    spike_forward = zeros(1,T);
    spikeL(round(1/dt)) = 1;
    spikeL_all_bout = spikeL;
    pos_bout = round(1/dt);
    s_bout = 1;
    
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
        
        
        if spikeL(i)>0 & rr(i)>rl(i);
            C(i+1) = C(i);
            Light_cur_R(i+1) = Itot*(1+C(i+1))/2;
            Light_cur_L(i+1) = Itot*(1-C(i+1))/2;
        end
        
        if spikeL(i)>0 & rl(i)>rr(i);
            C(i+1) = C(i);
            Light_cur_R(i+1) = Itot*(1+C(i+1))/2;
            Light_cur_L(i+1) = Itot*(1-C(i+1))/2;
        end
        if spikeL(i)==0;
            C(i+1) = C(i);
            Light_cur_R(i+1) = Itot*(1+C(i+1))/2;
            Light_cur_L(i+1) = Itot*(1-C(i+1))/2;
        end
        
        if rl(i)==rr(i);
            C(i+1) = C(i);
            Light_cur_R(i+1) = Itot*(1+C(i+1))/2;
            Light_cur_L(i+1) = Itot*(1-C(i+1))/2;
        end
        
    end
    
    dsig = rr-rl;
    
    for s_bout = 1:length(type_spikeL_all_bout)
        if type_spikeL_all_bout(s_bout)==0
            Sequence_delta_teta(s_bout) = sig_scout*randn(1);
        else
            Sequence_delta_teta(s_bout) = sig_turn*sign(dsig(time_bout(s_bout)))*abs(randn(1));
        end
    end
    
    % telegraph signal
    % transform into telegraph signal
    rlbin = rl;
    rrbin = rr;
    rlbin(rl>rr) = 1;
    rlbin(rl<rr) = 0;
    rrbin(rr>rl) = 1;
    rrbin(rr<rl) = 0;
    telsig = rlbin - rrbin;
    % find residence time
    ds = diff(telsig);
    stateonset_pos = find(ds>0); % loc of flips
    stateonset_neg = find(ds<0); % loc of flips
    if isempty(stateonset_neg)
        stateonset_neg=0;
    end
    if isempty(stateonset_pos)
        stateonset_neg=0;
    end
    restime_pos = [stateonset_pos(1) diff(stateonset_pos)]; % temps de résidence
    restime_neg = [stateonset_neg(1) diff(stateonset_neg)]; % temps de résidence
    restime_pos = restime_pos*dt; % en secondes
    restime_neg = restime_neg*dt; % en secondes
    restimex = [0:0.1:20];
    restimepospdf = computePDF(restimex, restime_pos);
    restimenegpdf = computePDF(restimex, restime_neg);
    
    %***
    subplot(2,1,1)
    hold on
    plot(restimex, smooth(restimepospdf,3))
    title('residence time of telegraph signal')
    ylabel('pdf')
    
    subplot(2,1,2)
    hold on
    plot(restimex, smooth(restimenegpdf,3))
    title('residence time of telegraph signal')
    ylabel('pdf')
    toc
    
    mean_bias = [mean_bias mean(Sequence_delta_teta)];
end

%%
figure
times=0:dt:dt*(T-dt);
hold on
plot(times,10000*spikeL,'r:','LineWidth',2)
plot(times,10000*spikeR,'b:','LineWidth',2)
plot(times,rl,'LineWidth',2,'color','r')
plot(times,rr,'LineWidth',2,'color','b')

% ---
% find residence time
ds = diff(telsig);
stateonset_pos = find(ds); % loc of flips
restime = [stateonset_pos(1) diff(stateonset_pos)]; % temps de résidence
restime = restime*dt; % en secondes 
restimex = [0:0.1:80];
restimepospdf = computePDF(restimex, restime);

%*** 
figure
plot(restimex, restimepospdf)
title('residence time of telegraph signal')
ylabel('pdf')

%***
figure
plot(c, mean_bias)
