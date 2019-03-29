function [boutindexes] = BoutSpotFromCoord(coordinates, framerate, checkplot)

% Extraction of bout indexes from x y coordinated and angle of the fish
% INPUT
% ---
% coordinates   : x and y coordinates of the mass center of the fish
%                 size(coordinates) = t, dim (=2), seq
%                 t : number of recorded time points
%                 dim : x, y
%                 seq : number of different sequences in experiment
% framerate     : framerate per sequence
% checkplot     : == 1 : plot trajectories and spotted bouts
%
% OUTPUT
% ---
% boutindexes   : indexes of bouts in the sequence

%--------------------------------------------------------------------------
coord = coordinates;
fHz = framerate;
coord(coord==0) = NaN;
    
x = squeeze(coord(:,1,:))';
y = squeeze(coord(:,2,:))';

% get the differentials and their squares
dx = diff(x, 1, 2);
dxcarr = dx.^2;

dy = diff(y, 1, 2);
dycarr = dy.^2;

% get variances
vardxy = nanvar(dx(:)+dy(:));

% get the significant displacement 
sigdisplacementmatrix = ((dxcarr'+dycarr')/vardxy)';
datasetSize = size(sigdisplacementmatrix,1);

%%
% tunable parameters 
minIPI = 0.5; %minimum inter-peak interval (in secs)
minh = prctile(sigdisplacementmatrix(:),80);% minimum peak height

boutwindow = 0.25;%secs
%% ---

% initialization loop to get the maximum number of peaks (bouts) per
% sequence

% ini
pkmaxnb = 0;
% /
for i = 1 : datasetSize
    pks = findpeaks(sigdisplacementmatrix(i,:), 'MinPeakHeight', minh,...
        'MinPeakDistance', minIPI*fHz(i));
    if pkmaxnb < length(pks)
        pkmaxnb = length(pks);
    end
end
%% ---

% ini
maxbout = NaN(size(sigdisplacementmatrix,1), pkmaxnb);
boutindexes = NaN(size(sigdisplacementmatrix,1),pkmaxnb);
% /
for i = 1 : datasetSize
    sdm = sigdisplacementmatrix(i,:);
    [pks, locs] = findpeaks(sdm, 'MinPeakHeight', minh,...
        'MinPeakDistance', minIPI*fHz(i));
    maxbout(i, 1:length(locs)) = locs;
    tau = floor(fHz(i)*boutwindow);
    if tau < 1
        tau = 1;
    end
    
    for j=1:length(locs)
        subreg = locs(j)-tau : locs(j)+tau;
        neg = length(find(subreg<=0));
        oversize = length(find(subreg>length(sdm)+1));
        if isempty(neg)
            neg = 0;
        end
        subreg = subreg(subreg > 0);
        if ~isempty(oversize)
            subreg = subreg(1:end-oversize);
        end
        prereg = locs(j)-tau+neg : locs(j)-1-floor(2*tau/3);
        
        subx = x(i, subreg); % subregion of the signal around the bout
        px = polyfit(prereg, x(i, prereg), 1); % linear fit of the pre-bout signal
        subreg_redressx = subx-px(2)-px(1)*subreg; % substract the resulting fitting function       
        sd = std(subreg_redressx(1:length(prereg)));
        binsubregx = subx;
        binsubregx( abs(subreg_redressx) > 3*sd) = 1;
        binsubregx( abs(subreg_redressx) <= 3*sd) = 0;

        suby = y(i, subreg);
        py = polyfit(prereg, y(i, prereg), 1); % linear fit of the pre-bout signal
        subreg_redressy = suby-py(2)-py(1)*subreg;
        subreg_redressy(1:3)=0;
        sd = std(subreg_redressy(1:length(prereg)));
        binsubregy = suby;
        binsubregy( abs(subreg_redressy) > 3*sd) = 1;
        binsubregy( abs(subreg_redressy) <= 3*sd) = 0;
        
        if nanstd(subx) < 0.5
            binsubregx(:) = 0;
        end
        if nanstd(suby) < 0.5
            binsubregy(:) = 0;
        end
        [~, transition] = find(diff([binsubregx; binsubregy],1, 2)==1,1); % first val beyond sig
        if ~isempty(transition)
            locs(j)=locs(j)+neg-tau+transition-1;
        else
            locs(j) = NaN;
        end
    end
    
    pks(isnan(locs)) = [];
    locs(isnan(locs)) = [];
    
    if checkplot == 1
        h = figure;
        subplot(1,2,1)
        plot(x(i,:))
        hold on
        plot(y(i,:))
        plot(sdm, 'k')
        plot(locs, x(i,locs),'o',locs, y(i,locs), 'o')
        subplot(1,2,2)
        plot3(1:length(x(i,:)),x(i,:),y(i,:))
        waitfor(h)
    end
    
    boutindexes(i,1:length(locs)) = locs;
end
