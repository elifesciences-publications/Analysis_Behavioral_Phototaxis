%% Simulation
%DD=[];
%for bias_slope=0:0.01:.1;
 
cor2=[];
w1=wturn;
w2=wfor;
p=pturn;
pflip=0.18;
bias_slope=1.5; % 

lrst=1;
dx=[];
x=0;
qq=[];
lr=[];
for t=1:100000
    xlast=wrapToPi(x(end));
    q=pflip*(1-0.5*lrst*bias_slope*stereo_bias(xlast));
    if q<0; q=0;
    end
    if q>1; q=1;
    end
    lrst=lrst*((rand>q)*2-1);
    turn=rand<pturn;
    if turn==1;
        dxt=lrst*abs(w1*randn);
    else dxt=w2*randn;
    end
    dx=[dx,dxt];
    x=[x,x(end)+dxt];
    qq=[qq, q];
    lr=[lr, lrst];
end
 
hist(wrapToPi(x),10)

 
% extract the drift
nbins=11
xwrap2pi=wrapToPi(x);
[Xbin,elts_per_bin, v2bin, ~, ~ ] = BinsWithEqualNbofElements(stereo_bias(xwrap2pi(1:end-1)), (dx),10, 12 );
Drift = mean(v2bin,2);
p = polyfit(Xbin,Drift,1)
Drift_fit= polyval(p,Xbin);
hold off
plot(Xbin,Drift, '--','Linewidth', 1.5);
hold on;
plot(Xbin, Drift_fit, '.-', 'LineWidth', 1.5)
%DD=[DD, p(1)]
%end
%pflip*bias_slope*sqrt(2/pi)*wturn*pturn/p(1)

 
hold off;

