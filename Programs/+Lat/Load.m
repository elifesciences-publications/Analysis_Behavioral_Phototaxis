
%% Load stereo-phototaxis data (lateralized experiments)

figpath = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Figures201808/';
path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
name = 'lateralized_exps.mat';
if ~exist([path name], 'file')
    prompt  = 'file non-existing. Create ? (y/n)';
    x = input(prompt, 's');
    if x == 'y'
        [El, DatesUsed, FishUsed] = Lat.lateralization_pool(path);
    else
        disp('canceled')
        return
    end
else
    load([path name], 'El')
    disp('El loaded')
end

%clearvars -except El DatesUsed FishUsed

threshold = 5;
[XLat,Xfilt,Xlab, FishID, Advancement, TimeBout, xCoord, yCoord, sequencesperfish, boutsperfish]...
    = Lat.delete_f_not_enough_bouts(El, threshold);

IBI = diff(TimeBout, 1, 2);
pxmm = 11.5;

% ----- Useful variables -----
% --- illumination ---
percPm =1;
lum_lin = luminosity_linear(percPm);

L = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),pi-abs(pi-mod(Xfilt,2*pi))); 
R = interp1(deg2rad(lum_lin(:,1)),lum_lin(:,2),abs(pi-mod(Xfilt,2*pi)));
DIlr = L-R;

% --- select the X on which to work on
X = Xfilt;
X(TimeBout>60) = NaN; % temporal criterium

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
dXl = diff(XLat, 1, 2 );
dXf  = diff(Xfilt, 1, 2 );

% --- per fish
different_fish = unique(FishID);
fishdXLstd = NaN(length(different_fish), 1);
fishdXLmean = NaN(length(different_fish), 1);
fishBiasL = NaN(length(FishID),1);
fishDevL = NaN(length(FishID),1);
scount = 1;
for i = different_fish'
    fishseqs = find(FishID == i);
    dx = dXl(fishseqs,:);
    fishdXLstd(i) = nanstd(dx(:));
    fishdXLmean(i) = nanmean(dx(:));
    
    %--- vectors same size as dataset ---
    fishBiasL(scount : scount+length(fishseqs)-1) = repmat(fishdXLmean(i), length(fishseqs), 1);
    fishDevL(scount : scount+length(fishseqs)-1) = repmat(fishdXLstd(i), length(fishseqs), 1);
    scount = scount + length(fishseqs);
end
dXLssbiais = dXl - repmat(fishBiasL, [1 size(dXl,2)]);
dXLn = dXLssbiais ./ repmat(fishDevL, [1 size(dXl,2)]);
dXLssdev = dXl ./ repmat(fishDevL, [1 size(dXl,2)]);
