%% Pool or Load pooled data from temporal phototaxis data
% on enucleated fish

pooled_data_path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
name = 'enucleated_exps.mat';
if ~exist([pooled_data_path name], 'file')
    prompt  = 'file non-existing. Create ? (y/n)';
    x = input(prompt, 's');
    if x == 'y'
        [Ee, DatesUsed, FishUsed] = Enuc.enuc_pool(pooled_data_path);
    else
        disp('canceled')
        return
    end
else
    load([pooled_data_path name], 'Ee')
    disp('Ee loaded')
end

%clearvars -except El DatesUsed FishUsed

XSource = Ee.AngleSource;
Xlab = Ee.AngleLab;
Xfilt = Ee.AngleSourceFiltered;
FishN = Ee.FishN;
Advancement = Ee.R;
Luminosity = Ee.Luminosity;
TimeBout = Ee.TimeBout;
T = Ee.T;
xCoord = Ee.xCoord;
yCoord = Ee.yCoord;
IBI = diff(TimeBout, 1, 2);
pxmm = 11.5;

% ----- Useful variables -----
% --- illumination ---
percPm =1;
lum_lin = luminosity_linear(percPm);

dLum = diff(Luminosity,1,2);

% --- select the X on which to work on
X = Xfilt;

% --- coordinates
dxCoord = diff(xCoord, 1, 2);
dyCoord = diff(yCoord, 1, 2);
Dist = sqrt(dxCoord.^2+ dyCoord.^2);

% --- advancement
trajOrientation = wrapToPi(atan(dyCoord./dxCoord));
trajOrientation(dxCoord < 0) = trajOrientation(dxCoord < 0) - pi;
Alpha = wrapToPi(Xlab(:, 1:end-1)) + trajOrientation;
relAdvancement = cos(Alpha);
Advancement = Dist.*cos(Alpha);

% --- dX
dX = diff(XSource, 1, 2 );

% --- correct wrong angles 
wa = find(Advancement<0 & abs(dX)<pi/2);
XSource(wa) =  NaN;
Xfilt(wa) = NaN;
dX  = diff(Xfilt, 1, 2 );
dX(abs(dX)>2*pi/3) = NaN;
%........................
% some stats on sequences
    % number of fish
seqPerFish = histc(FishN, unique(FishN));
NFish = length(seqPerFish);
meanNofSeqPerFish = mean(seqPerFish);

NofSeqs = size( XSource, 1 );
boutsPerSeq = size( XSource, 2 ) - sum(isnan(XSource),2);
maxnseq = max(seqPerFish);
