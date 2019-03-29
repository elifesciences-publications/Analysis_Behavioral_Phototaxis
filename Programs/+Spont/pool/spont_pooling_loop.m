function [Xi, TimeBouti, xCoordi, yCoordi, FishNi] = spont_pooling_loop(Dates, Fish, l, L, exptype, intmax, timemin, minBoutNumber)

%%
% --- for display ---
empty = cellfun(@isempty, Fish);
e = sum(empty(:));
nbOfExperiments = numel(Fish) - e;

Xi = NaN(l,L);
TimeBouti = NaN(l,L);
CoordinatesBouti = NaN(l,L,2);
FishNi = NaN(l, 1);

count = 1;
row = 1;
for p = 1: length(Dates)
    for q = 1: length(Fish(p,:)) 
        
        % --- select and load experiment
        if iscellstr(Fish(p,q)) == 0
            continue
        end
        disp(['experiment # ' num2str(count) '/' num2str(nbOfExperiments)])
        
        fish = char(Fish(p,q));        
        date = char(Dates(p));
        
        if intmax ~= 0
            datapath = ['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype filesep intmax filesep date filesep fish 'ini.mat'];
        else
            datapath = ['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype filesep date filesep fish 'ini.mat'];
        end
        if ~exist(datapath, 'file')
            continue
        end
        load(datapath, 'Dini');        
        
        % --- Remove the short sequences ---
        [angleLab_ini, coordinates_ini, framerateInfo_ini] = remove_sequence_ini(Dini, timemin);
        framerate_ini = framerateInfo_ini(:,3);
        
        % --- Get the cumulative angle ---
        %     (lab reference, the angle is initially wrapped to 360) 
        [angleLabCum_ini] = angle_cum(angleLab_ini, framerateInfo_ini);
        
        if isempty(angleLabCum_ini)
            continue
        end
        
        % --- Correct for potential aberrant bouts (head/tail misplacement)
        angleLabCumCorr_ini = find_aberrant_bouts(angleLabCum_ini, coordinates_ini, 0);
        
        % --- Get a smoothed version of the angle ---
        angleLabCumCorrSm_ini = smoothCoord(angleLabCumCorr_ini);
        
        % --- Spot real bouts
        fig = 0;
        boutIndices_ini = BoutSpot(coordinates_ini, angleLabCumCorr_ini, framerateInfo_ini, fig);

        % --- Get the angle before each bout ---
        [angleLabBout_ini, coordinatesBouti] = findAngBout(boutIndices_ini, angleLabCumCorrSm_ini, coordinates_ini);
        
        % --- Delete sequences with too few bouts ---
        BoutNumberPerSeq = size(angleLabBout_ini, 2) - sum(isnan(angleLabBout_ini), 2);
        seqToDel = find( BoutNumberPerSeq < minBoutNumber );
        if ~isempty(seqToDel)
            angleLabBout_ini(seqToDel, :) = [];
            coordinatesBouti(seqToDel, :, :) = [];
            framerate_ini(seqToDel) = [];
            boutIndices_ini(seqToDel,:) = [];
        end
        timeofbouts_ini = boutIndices_ini./framerate_ini;

        %store cum angles & lum
        %row = find( isnan(Xi(:,1)), 1, 'first');
        sizeExp = size(angleLabBout_ini);
        
        Xi( row : row + sizeExp(1)-1 , 1 : sizeExp(2)) = deg2rad(angleLabBout_ini);
        TimeBouti( row : row + sizeExp(1)-1 , 1 : sizeExp(2) ) = timeofbouts_ini;
        CoordinatesBouti( row : row + sizeExp(1)-1 , 1 : sizeExp(2), :) = coordinatesBouti;
        FishNi(row : row + sizeExp(1)-1) = ones(sizeExp(1), 1)*count;
       
        row = row +  sizeExp(1);
        disp(['row : ' num2str(row) ' size :' num2str(sizeExp(1))])
        
        count = count + 1;
    end
end

[xCoordi] = deleteEmptyRows(CoordinatesBouti(:,:,1));
[yCoordi] = deleteEmptyRows(CoordinatesBouti(:,:,2));
[Xi] = deleteEmptyRows(Xi);
[TimeBouti] = deleteEmptyRows(TimeBouti);
FishNi = deleteEmptyRows(FishNi);

%% delete sequences where the fish was head/tail misplaced

dxCoordi = diff(xCoordi, 1, 2);
dyCoordi = diff(yCoordi, 1, 2);
Dist = sqrt(dxCoordi.^2+ dyCoordi.^2);

% calculate advancement
trajOrientation = wrapToPi(atan(dyCoordi./dxCoordi));
trajOrientation(dxCoordi < 0) = trajOrientation(dxCoordi < 0) - pi;
dAlpha = wrapToPi(Xi(:, 1:end-1)) + trajOrientation;

Ri = Dist.*cos(dAlpha);

% completely delete sequences where cumulative advancement <0
cumulativeAdvancement = nansum(Ri, 2);
cumathresh = 100;
seqs_to_delete = find(cumulativeAdvancement < cumathresh);

Xi(seqs_to_delete, :) = []; 
TimeBouti(seqs_to_delete, :) = []; 
xCoordi(seqs_to_delete, :) = []; 
yCoordi(seqs_to_delete, :) = []; 
FishNi(seqs_to_delete, :) = []; 

%% replace aberrant points by NaNs
dXi = diff( Xi, 1, 2 );

dxCoordi = diff(xCoordi, 1, 2);
dyCoordi = diff(yCoordi, 1, 2);
Dist = sqrt(dxCoordi.^2+ dyCoordi.^2);

% calculate advancement
trajOrientation = wrapToPi(atan(dyCoordi./dxCoordi));
trajOrientation(dxCoordi < 0) = trajOrientation(dxCoordi < 0) - pi;
dAlpha = wrapToPi(Xi(:, 1:end-1)) + trajOrientation;

Ri = Dist.*cos(dAlpha);

% delete aberrant points
[S, B] = find(Ri<0 & abs(dXi)<pi/2);
Xi(S,B) = NaN;
TimeBouti(S,B) = NaN;
xCoordi(S,B) = NaN;
yCoordi(S,B) = NaN;

%% finally clean up all matrices
% delete empty rows
[Xi, delrows] = deleteEmptyRows(Xi);
[xCoordi] = deleteEmptyRows(xCoordi);
[yCoordi] = deleteEmptyRows(yCoordi);
TimeBouti(delrows,:) = [];

FishNi(delrows) = [];
