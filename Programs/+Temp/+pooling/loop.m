function [Xrad, Xlab, Lu, TimeBout, xCoord, yCoord, FishN] = loop(Dates, Fish, l, L, exptype, intmax, timemin, minBoutNumber)
% Loop for pooling temporal phototaxis experiments
% 2018/12/04

%%
% --- for display ---
empty = cellfun(@isempty, Fish);
e = sum(empty(:));
nbOfExperiments = numel(Fish) - e;

Xrad = NaN(l,L);
Xlab = NaN(l,L);
Lu = NaN(l,L);
TimeBout = NaN(l,L);
CoordinatesBout = NaN(l,L,2);
FishN = NaN(l, 1);

count = 1;
for p = 1: length(Dates)
    for q = 1: length(Fish(p,:)) 
        
        % --- select and load experiment
        if iscellstr(Fish(p,q)) == 0
            continue
        end
        disp(['experiment # ' num2str(count) '/' num2str(nbOfExperiments)])
        
        fish = char(Fish(p,q));        
        date = char(Dates(p));
        
        load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype '/' intmax '/' date '/' fish '/data' fish '.mat']);
        
        % ----- Remove the short sequences
        [angleLab ,angleSource, ~,framerateInfo ,coordinates ,luminosity] = remove_sequence(D,timemin);
        framerate = framerateInfo(:,3);        

        % --- Correct for potential aberrant bouts (head/tail misplacement)
        [angleLabCum] = angle_cum(angleLab, framerateInfo);
        angleSourceCorr = find_aberrant_bouts(angleSource, coordinates, 0);
        
        % --- Spot real bouts
        fig = 0;
        boutIndices = BoutSpot(coordinates, angleSourceCorr, framerateInfo, fig);
        [angleSourceCorrSm, ~, ~] = smoothCoord(angleSourceCorr);
        [angleLabCumSm, ~, ~] = smoothCoord(angleLabCum);
        
        % ----- and illumination, angle and coordinates before bouts
        [lumBout, angleBoutSource, angleBoutLab, coordinatesBout] = Temp.findAngleLumBout(boutIndices, luminosity, angleSourceCorrSm, angleLabCumSm, coordinates);
        timeofbouts = boutIndices./framerate;
        
        % --- delete sequences with too few bouts ---
        BoutNumberPerSeq = size(angleBoutSource, 2) - sum(isnan(angleBoutSource), 2);
        seqToDel = find( BoutNumberPerSeq < minBoutNumber );
        if ~isempty(seqToDel)
            angleBoutSource(seqToDel, :) = [];
            lumBout(seqToDel, :) = [];
            coordinatesBout(seqToDel, :, :) = [];
            timeofbouts(seqToDel, :) = [];
            angleBoutLab(seqToDel, :) = [];
        end

        %store cum angles & lum
        row = find(isnan(Xrad(:,1)),1, 'first');
        sizeExp = size(angleBoutSource);
        Xrad( row : row + sizeExp(1)-1 , 1 : sizeExp(2)) = deg2rad(angleBoutSource) ;
        Xlab( row : row + sizeExp(1)-1 , 1 : sizeExp(2)) = deg2rad(angleBoutLab);
        Lu( row : row + sizeExp(1)- 1 , 1 : sizeExp(2)) =  lumBout;
        TimeBout( row : row + sizeExp(1)-1 , 1 : sizeExp(2) ) = timeofbouts;
        CoordinatesBout( row : row + sizeExp(1)-1 , 1 : sizeExp(2) , :) = coordinatesBout;
         
        FishN(row : row + sizeExp(1)-1) = ones(sizeExp(1), 1)*count;
        
        count = count + 1;
    end
end

[xCoord] = (CoordinatesBout(:,:,1)); % deleteEmptyRows(CoordinatesBout(:,:,1));
[yCoord] = (CoordinatesBout(:,:,2)); %deleteEmptyRows

%% delete sequences where the fish was head/tail misplaced

dxCoord = diff(xCoord, 1, 2);
dyCoord = diff(yCoord, 1, 2);
Dist = sqrt(dxCoord.^2+ dyCoord.^2);

% calculate advancement
trajOrientation = wrapToPi(atan(dyCoord./dxCoord));
trajOrientation(dxCoord < 0) = trajOrientation(dxCoord < 0) - pi;
dAlpha = wrapToPi(Xlab(:, 1:end-1)) + trajOrientation;

R = Dist.*cos(dAlpha);

% completely delete sequences where cumulative advancement <0
cumulativeAdvancement = nansum(R, 2);
cumathresh = 100;
seqs_to_delete = find(cumulativeAdvancement < cumathresh);

Xrad(seqs_to_delete, :) = []; 
Xlab(seqs_to_delete, :) = []; 
Lu(seqs_to_delete, :) = []; 
TimeBout(seqs_to_delete, :) = []; 
xCoord(seqs_to_delete, :) = []; 
yCoord(seqs_to_delete, :) = []; 
FishN(seqs_to_delete, :) = []; 

%% replace aberrant points by NaNs
dX = diff( Xrad, 1, 2 );

dxCoord = diff(xCoord, 1, 2);
dyCoord = diff(yCoord, 1, 2);
Dist = sqrt(dxCoord.^2+ dyCoord.^2);

% calculate advancement
trajOrientation = wrapToPi(atan(dyCoord./dxCoord));
trajOrientation(dxCoord < 0) = trajOrientation(dxCoord < 0) - pi;
dAlpha = wrapToPi(Xlab(:, 1:end-1)) + trajOrientation;

R = Dist.*cos(dAlpha);

% delete aberrant points
[S, B] = find(R<0 & abs(dX)<pi/2);
Xrad(S,B) = NaN;
Xlab(S,B) = NaN;
Lu(S,B) = NaN;
TimeBout(S,B) = NaN;
xCoord(S,B) = NaN;
yCoord(S,B) = NaN;

%% finally clean up all matrices
% delete empty rows
[Xrad, delrows] = deleteEmptyRows(Xrad);
[Xlab] = deleteEmptyRows(Xlab);
[Lu] = deleteEmptyRows(Lu);
[xCoord] = deleteEmptyRows(xCoord);
[yCoord] = deleteEmptyRows(yCoord);
[TimeBout] = deleteEmptyRows(TimeBout);

FishN(delrows) = [];
