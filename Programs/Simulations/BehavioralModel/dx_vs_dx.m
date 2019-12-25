function[x, C1, convC1err, t, Cor, MSD] = dx_vs_dx(wturn, wfor, pturn, p_flip)

% compute dX(n) as a function of dX(n-1)
% consider two Markov chain for forward to turn and Right to left

% values from the turn distribution
%pturn=0.42;        % ratio of turning bouts 0.465
%wturn=0.59;        % std of the turn distribution in radian
%wfor=0.092;        % sdt of the forward distribution in radian
%p_flip= 0.1889;    % probability of flipping direction (R->L or L->R) - to be fitted


afor=1-pturn;   % ratio of forward bouts
p_TF = 1-pturn;   % probability of going from turn to forward state

                       
% creates the distributions
x=[-2:0.01:2];
f = afor*normpdf(x,0,wfor);
g = pturn*normpdf(x,0,wturn);
h = g./(f+g); % probability of being in turn state
alpha = sum(g)/sum(f); % ratio of turn vs forward bouts
                       % p_FT = alpha*p_TF


%% get error estimation

[errfit] = theta_measure_error_estimation(x);

%% calculate <dXn> as a function of <dXn-1>
C1=sqrt(2/pi)*(1-2*p_flip)*sign(x).*(h*(1-p_TF)*wturn);

%--- convolve with the gaussian kernel from error fit ---
convC1err = conv(C1,errfit, 'same')/sum(errfit);

%***
%figure;
plot(x,smooth(C1));
hold on
plot(x,convC1err);

%% calculate autocorrelation function
t=[0:30];
Cor=2/pi*(1-2*p_flip).^t*((1-p_TF)*wturn)^2/((1-p_TF)*wturn^2+p_TF*wfor^2);
Cor(1)=1;
figure;
plot(t,Cor)

%% calculate MSD
[d1,d2]=meshgrid(t,t);
MM=zeros(size(t,2), size(t,2));
for p1=1:size(t,2)
    for p2=1:size(t,2)
        MM(p1,p2)=Cor((abs(p1-p2)+1));
    end
end
MSD=[0];
for p=1:(size(t,2)-1)
    MMpart=MM(1:p,1:p); 
    MSD=[MSD sum(MMpart(:))];
end
MSD=((1-p_TF)*wturn^2+p_TF*wfor^2)*MSD;
%figure;
plot(t,MSD);

