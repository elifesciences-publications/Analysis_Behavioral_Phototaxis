% get IBI exp;
hold off

 
load('/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/spontaneous_swim.mat');
texp=Es.TimeBout;
inter_bout=texp(:,2:end)-texp(:,1:end-1);
inter_bout=inter_bout(:);
w=find(isfinite(inter_bout) & inter_bout<5); % get rid of very long IBI
IBI_exp=inter_bout(w);

 
binsize=.4; % bin size is 333ms 
[IBI_hist,edges]=histcounts(IBI_exp,[0:binsize:20]);
binval=edges(1:end-1)+binsize/2;
plot(binval,IBI_hist/sum(IBI_hist)/binsize,'-*');

rflip = binsize*0.19/median(IBI_exp)
 
%% Numerical generation of bouts
prep=IBI_hist(2)/length(IBI_exp);
pbout=IBI_hist(3)/length(IBI_exp)/(1-prep);
% this value is to force the second and third point to be similar
trans1 = [prep,1-prep;pbout,1-pbout];
emis1=[0,1;1,0];

 
N=10000000;
[seq1,states1] = hmmgenerate(N,trans1,emis1);
%plot(states(1:50));
tnum=find(seq1==2);
IBI_num=tnum(:,2:end)-tnum(:,1:end-1);
[IBInum_hist,~]=histcounts(IBI_num,edges/binsize);%[0:1:60]);

 
figure;
subplot(2,1,1)
hold off; plot(binval,IBI_hist/sum(IBI_hist)/binsize,'-*', 'DisplayName', 'exp');
%[IBI_hist,~]=histcounts(IBI_exp,[0:0.1:20]);
%hold off; plot([0.05:0.1:19.95],IBI_hist/sum(IBI_hist)/0.1,'-*');
hold on; plot(binval,IBInum_hist/sum(IBInum_hist)/binsize,'-*', 'DisplayName', 'num');
xlabel('time in seconds');
ylabel('PDF(IBI)');
xlim([0,5]);
legend
 
subplot(2,1,2)
hold off; plot(binval,IBI_hist/length(IBI_exp),'-*', 'DisplayName', 'exp');
hold on; plot(binval,IBInum_hist/length(IBI_num),'-*', 'DisplayName', 'num');
xlabel('time in seconds');
ylabel('PDF(IBI)');
set(gca, 'YScale', 'log')
xlim([0 5])
legend
 
mean(IBI_num(find(IBI_num<15)))*binsize
mean(IBI_exp(find(IBI_exp<5)))
% we do not get exact same mean values, but the fit is correct

 
% Numerical generation of left and right states
%rflip=0.064; %0.09 adjust such that correlation match data
trans2 = [1-rflip,rflip;rflip,1-rflip];
emis2=[0,1; 1,0];
[seq2,states2] = hmmgenerate(N,trans2,emis2);

 
% Numerical generation of scoot and turn states
pturn=0.41;        % fraction of turning events
pt2f=0.47;         % this value is set by dx2n-1 vs dx2n fitting
pf2t=pt2f*pturn/(1-pturn);
pt2f=1-pturn;
pf2t=pturn;
trans3=[1-pt2f,pt2f;pf2t,1-pf2t];
emis3=[0,1; 1,0];
[seq3,states3] = hmmgenerate(N,trans3,emis3);

 
% generate angular sequence
wturn=0.59;  % std of the turn distribution in radian
wfor=0.092;   % sdt of the forward distribution in radian

 
seq=[seq1;seq2;seq3];
Tbouts=find(seq1==2);
%seq_bouts=seq(:,Tbouts);
ind_turn=find(seq(3,:)==2 & seq(1,:)==2);
ind_for=find(seq(3,:)==1 & seq(1,:)==2);
angle_turns=randn(1,length(ind_turn))*wturn;
angle_for=randn(1,length(ind_for))*wfor;

 
angles=zeros(N,1);
angles(ind_for)=angle_for;
sign_angle_turn=[seq2(ind_turn)>1.5]*2-1;
angles(ind_turn)=sign_angle_turn.*abs(angle_turns);
angles_disc=angles(Tbouts);

 
figure;
a=autocorr(angles_disc, 20);
plot(a(2:8),'*-')
hold on
plot(ACinormpf(2:8), 'o-')

xlabel('bout number');
ylabel('angles correlation');

 
% normalized correlation of thresholded angles as a function of IBI
thresh=0.22;
%angles_turn=angles(ind_turn);
ind_thresh=find(abs(angles)>thresh);
%ind_thresh=ind_turn;
angles_thresh=angles(ind_thresh);
IBI_turn=ind_thresh(2:end)-ind_thresh(1:end-1);
corr=angles_thresh(2:end).*angles_thresh(1:end-1)./abs(angles_thresh(2:end))./abs(angles_thresh(1:end-1));
[X,Y]=histprofile(IBI_turn,corr,50);

 
figure;
plot(X*binsize,Y,'*'); xlim([0,10]);
hold on

xlabel('IBI in seconds');
ylabel('sign(theta_n+1) * sign(theta_n)');
plot(X*binsize,exp(-X*2*rflip)); % this is the theoretical fit 
% should be perfect if only turn bouts were selected

%% compare to data and fit
turn_thresh = 0.22;
dX = dXissbiais;

[fig1, fig2, binvals1, mv1, stdv1] = Spont.autocorrelationVSinterboutinterval(dX, IBIi, turn_thresh, 'data ' );
close(fig2)

% simu
vq = interp1(X*binsize,Y,binvals1);

%theoretical fit 
theox = [0:0.01:binvals1(end)/binsize]*binsize;
theoy = exp(-[0:0.01:binvals1(end)/binsize]*2*rflip);

hold on
%plot(binvals1,vq,'*', 'DisplayName', 'simu'); 
plot(theox,theoy, 'Linewidth', 1.5, 'DisplayName', 'theoretical fit'); % this is the theoretical fit 
% should be perfect if only turn bouts were selected

xlim([0,binvals1(end)]);
xlabel('IBI in seconds');
ylabel('sign(theta_n+1) * sign(theta_n)');
