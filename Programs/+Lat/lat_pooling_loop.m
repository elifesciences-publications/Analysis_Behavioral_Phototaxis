function [XLat, Xlab, Xfilt, TimeBout, xCoord, yCoord, FishN, DatesUsed, FishUsed, CorrespondingFishN] ...
    = lat_pooling_loop(Dates, Fish, l, L, exptype, timemin, minBoutNumber)

XLat = NaN(l,L);
Xlab =  NaN(l,L);
Xfilt = NaN(l,L);
TimeBout = NaN(l,L);
CoordinatesBout = NaN(l,L,2);
FishN = NaN(l, 1);

DatesUsed = cell(size(Dates));
FishUsed = cell(size(Fish));
CorrespondingFishN = NaN(size(Fish));

%%
count = 1;
for p = 1: length(Dates) 
    for q = 1: length(Fish(p,:)) 
        
        % --- load experiment
        if iscellstr(Fish(p,q)) == 0
            continue
        end
        fish = char(Fish(p,q));        
        date = char(Dates(p));
        disp([date fish])
        
        datapath = ['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype filesep date filesep fish '.mat'];
        if ~exist(datapath, 'file')
            continue
        end
        load(datapath, 'E');
        
        % ----- Selection : on sequences length (short sequences removed)
        [angleLab, angleSource, angleLabFilt , framerateInfo, coordinates, ~] = remove_sequence(E, timemin);
        if isempty(angleSource)
            continue
        end
        framerate = framerateInfo(:,3);
        
        % --- Correction : for potential aberrant bouts (head/tail misplacement)
        fig = 0;
%         angleSourceCorr = find_aberrant_bouts(angleSource, coordinates, fig);
%         
%         % --- Correct angle relative to source
%         if ~isfield(E.experiment, 'Illum')
%             angleSourceCorr2 = angleSourceCorr + (angleLab(:,1) - angleSource(:,1));
%             angleLabFilt = angleLabFilt + (angleLab(:,1) - angleSource(:,1));
%         else
%             angleSourceCorr2 = angleSourceCorr;
%             angleLabFilt = smoothCoord(angleSourceCorr2);
%         end
        
        angleSourceCorr2 = angleSource;
        
        % --- Spot real bouts
        fig = 0;
        [boutIndices, ~, ~] = BoutSpot(coordinates, angleSourceCorr2, framerateInfo, fig);
        [angleSourceCorrSm, ~, ~] = smoothCoord(angleSourceCorr2);
        
        % ----- and illumination and angle before bouts
        [angleBoutSource] = findVarBout(boutIndices, angleSourceCorrSm, 2);
        [angleBoutLab] = findVarBout(boutIndices, angleLab, 2);
        [angleBoutLabFilt] = findVarBout(boutIndices, angleLabFilt, 2);
        [coordinatesBout] = findVarBout(boutIndices, coordinates, 3);
        timeofbouts = boutIndices./framerate;
        
        % --- delete sequences with too few bouts ---
        BoutNumberPerSeq = size(angleBoutSource, 2) - sum(isnan(angleBoutSource), 2);
        seqToDel = find( BoutNumberPerSeq < minBoutNumber );
        if ~isempty(seqToDel)
            angleBoutSource(seqToDel, :) = [];
            angleBoutLab(seqToDel, :) = [];      
            angleBoutLabFilt(seqToDel, :) = [];   
            coordinatesBout(seqToDel, :, :) = [];
            timeofbouts(seqToDel, :) = [];
        end
        
        %store cum angles
        sizeExp = size(angleBoutSource);
        [row, ~] = find(isnan(XLat),1, 'first');
        
        XLat( row : row + sizeExp(1)-1 , 1 : sizeExp(2)) = deg2rad(angleBoutSource);
        Xlab( row : row + sizeExp(1)-1 , 1 : sizeExp(2)) = deg2rad(angleBoutLab);
        Xfilt( row : row + sizeExp(1)-1 , 1 : sizeExp(2)) = deg2rad(angleBoutLabFilt);
        TimeBout( row : row + sizeExp(1)-1 , 1 : sizeExp(2) ) = timeofbouts;
        CoordinatesBout( row : row + sizeExp(1)-1 , 1 : sizeExp(2) , :) = coordinatesBout;
        
        FishN(row : row + sizeExp(1)-1) = ones(sizeExp(1), 1)*count;
        
        disp(row)
        
       DatesUsed(p) = Dates(p);
       FishUsed(p,q) = Fish(p,q);
       CorrespondingFishN(p,q) = count;
       count = count + 1;
    end
end

%% delete sequences where the fish was head/tail misplaced

[xCoord] = CoordinatesBout(:,:,1);
[yCoord] = CoordinatesBout(:,:,2);

dxCoord = diff(xCoord, 1, 2);
dyCoord = diff(yCoord, 1, 2);
Dist = sqrt(dxCoord.^2+ dyCoord.^2);

% calculate advancement
trajOrientation = wrapToPi(atan(dyCoord./dxCoord));
trajOrientation(dxCoord < 0) = trajOrientation(dxCoord < 0) - pi;
dAlpha = wrapToPi(Xlab(:, 1:end-1)) + trajOrientation;

Advancement = Dist.*cos(dAlpha);

%visualizeTrajectoryAndAdv([1:72]', Xlab, xCoord, yCoord, Dist, Advancement, diff(XLat,1,2), trajOrientation, 0)


% completely delete sequences where cumulative advancement <0
cumulativeAdvancement = nansum(Advancement, 2);
cumathresh = 1;
seqs_to_delete = find(cumulativeAdvancement < cumathresh);

XLat(seqs_to_delete, :) = []; 
Xlab(seqs_to_delete, :) = []; 
Xfilt(seqs_to_delete, :) = []; 
TimeBout(seqs_to_delete, :) = []; 
xCoord(seqs_to_delete, :) = []; 
yCoord(seqs_to_delete, :) = []; 
FishN(seqs_to_delete, :) = []; 

%% replace aberrant points by NaNs
dxCoord = diff(xCoord, 1, 2);
dyCoord = diff(yCoord, 1, 2);
Dist = sqrt(dxCoord.^2+ dyCoord.^2);

% calculate advancement
trajOrientation = wrapToPi(atan(dyCoord./dxCoord));
trajOrientation(dxCoord < 0) = trajOrientation(dxCoord < 0) - pi;
dAlpha = wrapToPi(Xlab(:, 1:end-1)) + trajOrientation;

Advancement = Dist.*cos(dAlpha);

% delete aberrant points
% [S, B] = find(Advancement<0 & abs(dAlpha)<pi);
% XLat(S,B) = NaN;
% Xfilt(S,B) = NaN;
% Xlab(S,B) = NaN;
% TimeBout(S,B) = NaN;
% xCoord(S,B) = NaN;
% yCoord(S,B) = NaN;

%% finally clean up all matrices
% delete empty rows
[XLat, delrows] = deleteEmptyRows(XLat);
[Xlab] = deleteEmptyRows(Xlab);
[Xfilt] = deleteEmptyRows(Xfilt);
[xCoord] = deleteEmptyRows(xCoord);
[yCoord] = deleteEmptyRows(yCoord);
[TimeBout] = deleteEmptyRows(TimeBout);

FishN(delrows) = [];
