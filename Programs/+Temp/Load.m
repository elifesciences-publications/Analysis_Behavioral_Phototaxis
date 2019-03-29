% Load temporal phototaxis data
%%
% --- load E ---
path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
file = 'temporal_exps';
load([path file], 'E')

% --- first ini loop ---
R = 0;
C = 0;
for e = 1 : size(E,2)
    [r, c] = size(E(e).AngleSource);
    disp([num2str(r) ' ' num2str(c)])
    R = R + r;
    C = max(C, c);
end

% --- ini vars ---
X = NaN(R, C);
Indiv = NaN(R, 1);
Lum = NaN(R, C);
Exp = NaN(R, 1);
Adv = NaN(R, C-1);
T = NaN(R, C-1);
TB = NaN(R, C);
D = NaN(R,C-1);

rp = 1;
for e = 1 : size(E,2)
   [r, c] = size(E(e).AngleSource);
   X(rp : rp + r -1, 1 : c ) = E(e).AngleSource;
   Indiv(rp : rp + r -1 ) = E(e).FishN;
   Lum(rp : rp + r -1, 1 : c ) = E(e).Lum;
   Adv(rp : rp + r -1, 1 : c-1 ) = E(e).R;
   T(rp : rp + r -1, 1 : c-1 ) = E(e).T;
   TB (rp : rp + r -1, 1 : c ) = E(e).TimeBout;
   Exp(rp : rp + r -1) = ones(r, 1)*e;
   D(rp : rp + r -1, 1 : c-1 ) = E(e).Dist;
   rp = rp+r;
end

dLum = diff(Lum, 1,2);
Xcustom = X;
%Xcustom(TB < 100) = NaN;
dX = diff(Xcustom,1,2);
IBI = diff(TB, 1, 2);
IBI(IBI<0) = NaN;

pxmm = 11.5;

%% bias per fish
f = abs(diff(Indiv));
f(f>0) = 1;
FishID = 1+[0 ; cumsum(f)]; 

bias = NaN(1, length(unique(FishID)));
Bias = [];
Stdev = [];
for i = 1 : length(unique(FishID))
    dX1f = dX(FishID==i,:);
    bias(i) = nanmean(dX1f(:));
    Bias= [Bias; ones(sum(FishID==i),1)*bias(i)];
    Stdev = [Stdev; ones(sum(FishID==i),1)*nanstd(dX1f(:))];
end        
dX_corr = dX - Bias;
dX_norm = dX_corr./Stdev;

