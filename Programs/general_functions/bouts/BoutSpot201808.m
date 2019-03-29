function [boutindexes, minh] = BoutSpot(coordinates, angle, framerate, checkplot)

% Extraction of bout indexes from x y coordinated and angle of the fish
% INPUT
% ---
% coordinates   : x and y coordinated of the mass center of the fish
% angle_source  : cumulative angle of the fish relative to the light source (in degrees)
% framerate     : framerate per sequence
% checkplot     : 0 or 1 whether you want a plot or not
%
% OUTPUT
% ---
% boutindexes   : indexes of bouts in the sequence

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
%minh = nanmedian(logsigdisplacementmatrix,2) + 2*nanstd(logsigdisplacementmatrix, 1, 2);
%minh(minh<0) = 0.1;
%minh(minh>1) = 0.5;
minh = ones(datasetSize,1)*0.1;
% [~, sequence_length] = max(isnan(logsigdisplacementmatrix),[], 2);
% coeff_function = @(x) max(85, min(95, 77+0.03*x));
% percentile = coeff_function(sequence_length);
% minh = prctile(logsigdisplacementmatrix, percentile, 2);
% minh = diag(minh);
 pkmaxnb = 0;
for i = 1 : datasetSize
    lsdm = logsigdisplacementmatrix(i,:);
%     if find(isnan(sdm),1) < 20*fHz(i) %6 seconds = 200 points
%         minh(i) = min(minh(i), 2);
%         minh(i) = max(minh(i), 0.2);
%     elseif find(isnan(sdm),1) > 20*fHz(i) || isempty(find(isnan(sdm),1)) 
%         minh(i) = min(minh(i), 3);
%         minh(i) = max(minh(i), 1.8);
%     end
%     minh(i) = max(minh(i), 0.01);
    pks = findpeaks(lsdm , 'MinPeakHeight', minh(i),...
        'MinPeakDistance', minIPI*fHz(i));
    if pkmaxnb < length(pks)
        pkmaxnb = length(pks);
    end
end
%%
maxbout = NaN(size(sigdisplacementmatrix,1),pkmaxnb);
boutwindow = 0.25;%secs
boutindexes = NaN(size(sigdisplacementmatrix,1),pkmaxnb);

for i = 1 : datasetSize
    lsdm = logsigdisplacementmatrix(i,:);
    [pks, locs] = findpeaks(lsdm, 'MinPeakHeight', minh(i),...
        'MinPeakDistance', minIPI*fHz(i));
    maxbout(i, 1:length(locs)) = locs;
    tau = floor(fHz(i)*boutwindow);
    
    for j = 1 : length(locs)
        subreg = locs(j) - tau : locs(j) + tau;
        neg = length(find(subreg<=0));
        oversize = length(find(subreg>length(lsdm)+1));
        if isempty(neg)
            neg = 0;
        end
        subreg = subreg(subreg > 0);
        if ~isempty(oversize)
            subreg = subreg(1:end-oversize);
        end
        prereg = locs(j)-tau+neg : locs(j)-1-floor(2*tau/3);
        
        subx = mx(i, subreg); % subregion of the signal around the bout
        px = polyfit(prereg, mx(i, prereg), 1); % linear fit of the pre-bout signal
        subreg_redressx = subx-px(2)-px(1)*subreg; % substract the resulting fitting function       
        sd = std(subreg_redressx(1:length(prereg)));
        binsubregx = double( abs(subreg_redressx) > 3*sd);

        suby = my(i, subreg);
        py = polyfit(prereg, my(i, prereg), 1); % linear fit of the pre-bout signal
        subreg_redressy = suby-py(2)-py(1)*subreg;
        subreg_redressy(1:3)=0;
        sd = std(subreg_redressy(1:length(prereg)));
        binsubregy = double( abs(subreg_redressy) > 3*sd);
        
        if nanstd(diff(subx)) < 0.3
            binsubregx(:) = 0;
        end
        if nanstd(diff(suby)) < 0.3
            binsubregy(:) = 0;
        end
        [~, transition] = find(diff([binsubregx; binsubregy],1, 2)==1,1); % first val beyond sig
        if ~isempty(transition)
            locs(j) = locs(j)+neg-tau+transition-1;
        else
            locs(j) = NaN;
        end
    end
    
    pks(isnan(locs)) = [];
    locs(isnan(locs)) = [];
    
    if checkplot == 1
        h = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1)
        yyaxis left
        plot((1:length(angle(i,:)))/fHz(i), x(i,:), '-')
        hold on
        plot((1:length(angle(i,:)))/fHz(i), mx(i,:), '-')
        hold on
        plot((1:length(angle(i,:)))/fHz(i), y(i,:), '-')
        plot((1:length(angle(i,:)))/fHz(i), my(i,:), '-')
        plot( (1:length(angle(i,:)))/fHz(i), (angle(i,:)), '-', 'Color', [0.5 0.5 0.5], 'Linewidth', 1)
        plot( (1:length(angle(i,:)))/fHz(i), (mtheta(i,:)), '-k' , 'Linewidth', 1)
        plot(locs/fHz(i), (mtheta(i,locs)), 'sq', 'Color', [0 0 0], 'MarkerSize', 2, 'MarkerFaceColor', [0.4 0.4 0.9])
        plot(locs/fHz(i), x(i,locs),'o',locs/fHz(i), y(i,locs), 'o')
        yyaxis right
        lsdm(lsdm<0) = 0;
        plot((1:length(angle(i,1:end-1)))/fHz(i), lsdm)
        title(['seq : ' num2str(i) ', ' num2str(minh(i))])
        
        subplot(1,2,2);
        histogram(lsdm,100)
        waitfor(h)
    end
    boutindexes(i,1:length(locs)) = locs;
end
