% compute dX(n) as a function of dX(n-1)
% consider two Markov chain for forward to turn and Right to left

 
% define the two distribution for turn and forward bout
aturn=0.465;
wturn=0.61;
afor=1-aturn;
wfor=0.1036;
p_TF = 0.29;   % probability of going from turn to forward state
x=[-2:0.01:2];
f = afor*normpdf(x,0,wfor);
g = aturn*normpdf(x,0,wturn);
h = g./(f+g); % probability of being in turn state
alpha = sum(g)/sum(f); % ratio of turn vs forward bouts
                       % p_FT = alpha*p_TF
p_flip=0.28; % probability of flipping direction (R->L or L->R)
 
% <dXn> as a function of <dXn-1>
C1=sign(x).*h*(1-p_TF)*(1-2*p_flip)*sqrt(2/pi)*wturn;
figure; plot(x,smooth(C1,20));

 
% <abs(dXn)> as a function of <dXn-1>
C2= (h*p_TF+(1-h)*(1-alpha*p_TF))*sqrt(2/pi)*wfor+...
    (h*(1-p_TF)+(1-h)*alpha*p_TF)*sqrt(2/pi)*wturn;
figure;
plot(x,smooth(C2,20)); xlim([0 2]); 
hold on; plot([0 2], ...
    [(p_TF*wfor+(1-p_TF)*wturn)*sqrt(2/pi) (p_TF*wfor+(1-p_TF)*wturn)*sqrt(2/pi)]);
hold off