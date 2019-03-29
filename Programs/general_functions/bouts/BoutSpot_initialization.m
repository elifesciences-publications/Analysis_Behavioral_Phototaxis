function [nbouts] = BoutSpot_initialization(coordinates, angle_source, framerate)

% Extraction of bout indexes from x y coordinated and angle of the fish
% INPUT
% ---
% coordinates   : x and y coordinated of the mass center of the fish
% angle_source  : cumulative angle of the fish relative to the light source (in degrees)
% framerate     : framerate per sequence
%
% OUTPUT
% ---
% nbouts   : max number of bouts in a sequence

%--------------------------------------------------------------------------
coord = coordinates;
fHz = framerate(:,3);
coord(coord==0) = NaN;
    
x = squeeze(coord(:,1,:))';
y = squeeze(coord(:,2,:))';

% smooth all variables (x, y, theta)
[mtheta, mx, my] = smoothCoord(angle, x, y);

% get the differentials and their squares
dx = diff(mx, 1, 2);
dxcarr = dx.^2;

dy = diff(my, 1, 2);
dycarr = dy.^2;

dtheta = diff(angle, 1, 2);
dthetacarr = dtheta.^2;

% get variances
vardth = nanvar(dtheta(:));
vardxy = nanvar(dx(:)+dy(:));

% get the significant displacement 
sigdisplacementmatrix = ((dthetacarr'/vardth).*((dxcarr'+dycarr')/vardxy))';
logsigdisplacementmatrix = log(sigdisplacementmatrix);
logsigdisplacementmatrix(~isfinite(logsigdisplacementmatrix)) = NaN;
datasetSize = size(sigdisplacementmatrix,1);

%%
minIPI = 0.3; %minimum inter-peak interval (in secs)
[~, sequence_length] = max(isnan(sigdisplacementmatrix),[], 2);
coeff_function = @(x) max(80, min(90, 77+0.03*x));
percentile = coeff_function(sequence_length);
minh = prctile(sigdisplacementmatrix, percentile, 2);
minh = diag(minh);

pkmaxnb = 0;
for i = 1 : datasetSize
    sdm = sigdisplacementmatrix(i,:);
    if find(isnan(sdm),1) < 10*fHz(i) %6 seconds = 200 points
        minh(i) = min(minh(i), 2);
        minh(i) = max(minh(i), 0.1);
    elseif find(isnan(sdm),1) > 20*fHz(i) || isempty(find(isnan(sdm),1)) 
        minh(i) = min(minh(i), 5);
        minh(i) = max(minh(i), 1);
    end
    minh(i) = max(minh(i), 0.01);
    pks = findpeaks(sdm , 'MinPeakHeight', minh(i),...
        'MinPeakDistance', minIPI*fHz(i));
    if pkmaxnb < length(pks)
        pkmaxnb = length(pks);
    end
end

nbouts = pkmaxnb;