%% fit dx_n+1^2 as a function of dx_n^2
aturn=0.41; % ratio of turning bouts
wturn=0.59;  % std of the turn distribution in radian
wfor=0.092;   % sdt of the forward distribution in radian
wt=wturn^2;
wf=wfor^2;

% probability of flipping direction (R->L or L->R) - to be fitted
p_flip=0.19; 
 
afor=1-aturn;   % ratio of forward bouts
p_TF = 0.535;   %0.535 % probability of going from turn to forward state

% creates the distributions
x=[-2:0.01:2];
f = afor*normpdf(x,0,wfor);
g = aturn*normpdf(x,0,wturn);
h = g./(f+g);

 
alpha=aturn/(1-aturn);
k2=0.42;%0.45
vardt=aturn*wt+(1-aturn)*wf;

 
dt2=wf+k2*alpha*(wt-wf)+h*(wt-wf)*(1-k2-k2*alpha);
plot(x.^2,dt2/vardt, 'Linewidth', 1.5, 'DisplayName', 'fit');
xlim([0.001,2]);
ylim([0,2]);

%% fit data
autocorrXsq = nanmean(interp_mV1);

myfittype = fittype('amplitude_fit(x, wturn, wfor, aturn, k2)',...
            'problem', {'wturn', 'wfor', 'aturn'}, 'coefficients', {'k2'});
myfit = fit(allbinvals,autocorrXsq',myfittype,'StartPoint', [0.45],'problem', {wturn, wfor, aturn});

k2 = myfit.k2;
k2ci = confint(myfit,0.99);
dt2=(wf+k2*alpha*(wt-wf)+h*(wt-wf)*(1-k2-k2*alpha))/vardt;
dt2ci=(wf+k2ci(1)*alpha*(wt-wf)+h*(wt-wf)*(1-k2ci(1)-k2ci(1)*alpha))/vardt;

%***
% randomized mean
shadedErrorBar(allbinvals, nanmean(interp_mV1_rand), nanstd(interp_mV1_rand),...
    'lineprops',{'.','Color', [0.5 0.5 0.5], 'DisplayName', '<\delta\theta^2_{randomized}>'})
%real mean
shadedErrorBar(allbinvals, nanmean(interp_mV1), nanstd(interp_mV1),...
    'lineprops',{'-','Color', [0.2 0.4 0.9], 'LineWidth', 1.5, 'DisplayName', '<\delta\theta^2_{n+1}>'})
%simulation
shadedErrorBar(x.^2, dt2, -dt2+dt2ci,...
    'lineprops',{'-','Color', [0.9 0.4 0.2], 'LineWidth', 1.5, 'DisplayName', 'fit with .99 CI'})

ylabel('<\delta\theta^2_{n+1}>')
xlabel('<\delta\theta^2_{n}>')
legend

ax = gca;
ax.XScale = 'log';
ax.FontName = 'TimesNewRoman';
ax.FontSize = 14;