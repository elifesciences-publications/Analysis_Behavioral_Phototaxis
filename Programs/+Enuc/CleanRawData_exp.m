% clean D before analysis

%% Enuc experiment

[Dates, Fish, exptype] = Enuc.DatesFish();
timemin = 5;

%% load previously stored parameters
data_root = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/';
load([data_root 'EnucExperimentsUsed.mat'], 'EnucExperimentsUsedAll');
EnucExperimentsUsed = EnucExperimentsUsedAll;
clear EnucExperimentsUsedAll

%%

EnucExperimentsUsedNew = cell(numel(Fish),3);
count = 1;

for p = 1: length(Dates) 
    for q = 1: length(Fish(p,:)) 
        
        if iscellstr(Fish(p,q)) == 0
            continue
        end
        fish = char(Fish(p,q));        
        date = char(Dates(p));
        
        checkeddates = find(strcmp(date, EnucExperimentsUsed(:,1)));
        if ~isempty(checkeddates)
            checkedfish = find(strcmp(fish, EnucExperimentsUsed(checkeddates,2)));
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
        [angleLab, angleSource, angleLabFilt , framerateInfo, coordinates, ~] = remove_sequence(D, timemin);
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
        angleSourceCorr2 = angleSourceCorr;
        angleLabFilt = smoothCoord(angleSourceCorr2);
        
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
        angleLabFilt(del,:)=[];
        framerateInfo(del,:)=[];
        coordinates(:,:,del)=[];
        
        E = D;
        E.experiment.angle = angleLab;
        E.experiment.angleCum = angleSourceCorr2;
        E.experiment.angleFiltered = angleLabFilt;
        E.experiment.framerate = framerateInfo;
        E.experiment.coordinates = coordinates;
        
        EnucExperimentsUsedNew{count, 1} = date;
        EnucExperimentsUsedNew{count, 2} = fish;
        EnucExperimentsUsedNew{count, 3} = del;
        
        count = count +1;
        
        save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype '/' date '/' fish '.mat'], 'E');
        
    end
end

%%
cellcleaning = find(sum(cellfun(@isempty, EnucExperimentsUsedNew),2)==3);
EnucExperimentsUsedNew(cellcleaning,:) = [];
EnucExperimentsUsedAll = [EnucExperimentsUsed; EnucExperimentsUsedNew];
save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/' exptype '/EnucExperimentsUsed.mat'], 'EnucExperimentsUsedAll'); 
