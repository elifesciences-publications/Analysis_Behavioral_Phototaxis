% Generalized linear model regression

path = '/Users/karp/Documents/PhD/Reports/2018_BehavioralAnalysis/matfiles/';
file = 'lumangle';
load([path file], 'E')

e = 1;

X = E(e).AngleSource;
dX = diff(X,1,2);
Lu = E(e).Lum;
dLu = diff(Lu,1,2);

clear E path file

%% linearize responses and predictors

% responses : dXn
y = [dX, NaN(size(dX,1), 1)]';
y = y(:);

% for Xn and Lu : 
Xn = X'; 
Xn = Xn(:);

Lun = Lu'; 
Lun = Lun(:);

% add a first NaN column to dX and dLu
dXn_1 = [NaN(size(dX,1), 1), dX]';
dXn_1 = dXn_1(:);

dLu_1 = [NaN(size(dX,1), 1), dLu]';
dLu_1 = dLu_1(:);

%%

P = [Xn, dXn_1, Lun, dLu_1];

b = glmfit(P, y, 'normal', 'link', 'logit');

mdl = fitglm(P, y, 'linear', 'link', 'logit')
