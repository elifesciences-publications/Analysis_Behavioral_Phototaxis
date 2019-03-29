% clean D before analysis

%% Lat

Lat.DatesFish

%% load previously stored parameters
load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
    exptype '/LatExperimentsUsed.mat']); 
LatExperimentsUsedOld = LatExperimentsUsedAll;
clear LatExperimentsUsedAll
%%

LatExperimentsUsedNew = cell(numel(Fish),3);
count = 1;

for p = 1: length(Dates) 
    for q = 1: length(Fish(p,:)) 
        
        if iscellstr(Fish(p,q)) == 0
            continue
        end
        fish = char(Fish(p,q));        
        date = char(Dates(p));
        
        checkeddates = find(strcmp(date, LatExperimentsUsedOld(:,1)));
        if ~isempty(checkeddates)
            checkedfish = find(strcmp(fish, LatExperimentsUsedOld(checkeddates,2)));
            if ~isempty(checkedfish)
                continue
            end
        end
        
        load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype '/' date '/' fish '/data' fish '.mat'], 'D');
        disp([date ' ' fish])
        
        % ----- Selection : on framerate (experiment removed if too slow) ---
        if sum(D.experiment.framerate(:,3) < 20) > 1
            continue
        end
        
        % ----- Selection : on sequences length (short sequences removed)
        [angleLab, angleSource, angleFilt , framerateInfo, coordinates, ~] = remove_sequence(D, timemin);
        if isempty(angleSource)
            disp('to delete : ')
            disp([date fish])
            continue
        end
        framerate = framerateInfo(:,3);
        
        % --- Correction : for potential aberrant bouts (head/tail misplacement)
        fig = 0;
        angleSourceCorr = find_aberrant_bouts(angleSource, coordinates, fig);
        
        % --- Correct angle relative to source
        if ~isfield(D.experiment, 'Illum')
            angleSourceCorr2 = angleSourceCorr + (angleLab(:,1) - angleSource(:,1));
            angleFilt = angleFilt + (angleLab(:,1) - angleSource(:,1));
        else
            angleSourceCorr2 = angleSourceCorr;
            angleFilt = smoothCoord(angleSourceCorr2);
        end
        
        % --- Spot real bouts
        fig = 1;
        [boutIndices, ~, ~] = BoutSpot(coordinates, angleSourceCorr2, framerateInfo, fig);
        
        del = input('sequences to delete manually : ');
        
        if length(del) == size(angleSourceCorr,1)
            disp('whole experiment discarded : ')
            disp([fish date])
            continue
        end
        angleLab(del,:)=[];
        angleSourceCorr2(del,:)=[];
        angleFilt(del,:)=[];
        framerateInfo(del,:)=[];
        coordinates(:,:,del)=[];
        
        E = D;
        E.experiment.angle = angleLab;
        E.experiment.angleCum = angleSourceCorr2;
        E.experiment.angleFiltered = angleFilt;
        E.experiment.framerate = framerateInfo;
        E.experiment.coordinates = coordinates;
        
        LatExperimentsUsedNew{count, 1} = date;
        LatExperimentsUsedNew{count, 2} = fish;
        LatExperimentsUsedNew{count, 3} = del;
        
        count = count +1;
        
        save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype '/' date '/' fish '.mat'], 'E');
        
    end
end

%%
cellcleaning = find(sum(cellfun(@isempty, LatExperimentsUsedNew),2)==3);
LatExperimentsUsedNew(cellcleaning,:) = [];
LatExperimentsUsedAll = [LatExperimentsUsedOld; LatExperimentsUsedNew];
save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
    exptype '/LatExperimentsUsed.mat'], 'LatExperimentsUsedAll'); 
