function[DX1, DY1, t, ac, MSDsim] = simulate_markov(wturn, wfor, pturn, p_flip)

% Simulate spontaneous nagivation as a markov process

% wturn std of the turn bouts
% wfor std of the forward bouts
% pturn ratio of turning bouts
% p_flip probability of flipping direction (R->L or L->R) 

trans1 = [1-p_flip,p_flip;
         p_flip,1-p_flip];
% first state is left, second states is right
     
trans2 = [1-pturn,pturn;
          1-pturn,pturn];
% first state is forward, second states is turn


N=10000000;  % number of time steps
r = abs(normrnd(0,1,[1,N])); % normally distributed sequence

% generate states
emis = [0; 0];
[~,states1] = hmmgenerate(N,trans1,emis);
[~,states2] = hmmgenerate(N,trans2,emis);
states=[states1',states2'];


% index -1 for leftward, +1 for rightward
wL=find(states1 == 1);
wR=find(states1 == 2);
states1(wL)=-1; states1(wR)=1;

% states 2 gets the variance, with random sign if forward bout
wF=find(states2 == 1);
wT=find(states2 == 2);
states2(wF)=wfor; 
states2(wT)=wturn; 
states1(wF)=(rand([1,size(wF)])>0.5)*2-1;


% produce simulated reorientation and orientation sequence
dXsim=states2.*r.*states1;
Xsim=cumsum(dXsim);

%% calculate <dXn> as a function of <dXn-1>
nbins=100;
dXn2=dXsim(2:end);
dXn1=dXsim(1:end-1);
%[DX1,DY1, ~, ~ ] = histprofile((dXn1), (dXn2),nbins );
[DX1, elts_per_bin, v2bin, ~, ~] = BinsWithEqualNbofElements((dXn1), (dXn2),nbins, nbins+10 );
DY1 = mean(v2bin,2);
%hold off;

figure;
plot(DX1,DY1);
xlim([-1.5 1.5])
%% autocorr
t=[0:30];
AC = xcorr(dXsim);
ac = AC(end-size(dXsim,2)+1 : end-size(dXsim,2)+t(end)+1);
ac = ac/max(ac);

figure;
plot(t, ac);

%%
MSDsim=zeros(size(t,2),0);
for tt=t;
    val=var(Xsim((tt+1):end)-Xsim(1:end-tt));
    MSDsim(tt+1)=val;
end
figure;
plot(t,MSDsim);
    

