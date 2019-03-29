% ENUC
%% Experiments with enucleated samples :
%  clean raw data before analysis

%% Spont
timemin = 5;

%% load previously stored parameters
data_root = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/';
load([data_root 'EnucSpontUsed.mat'], 'EnucSpontUsedNew');
EnucSpontUsed = EnucSpontUsedAll;
clear EnucSpontUsedAll
%%
[Dates, Fish, exptype] = Enuc.DatesFish();

EnucSpontUsedNew = cell(50,3);
count = 1;
L = 0 ;
l = 0 ;

for p = 1: length(Dates)
    for q = 1: length(Fish(p,:))
        if iscellstr(Fish(p,q)) == 0 
            continue
        end
        fish = char(Fish(p,q));
        date = char(Dates(p));
        
        load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
                exptype '/' date '/' fish '/data' fish '.mat'], 'D');
        disp('loaded : ')
        disp([date ' ' fish])
        
        % check INITIALIZATION
        %__________________________________________________________________
        
        % --- 1st selection : sequences length (short sequences removed)
        % .............................................................
        [angle, coordinates, framerateInfo, del_length] = remove_sequence_ini(D, timemin);
        if ~isempty(del_length)
            if size(del_length,2)> size(del_length,1)
                del_length = del_length';
            end
        end
        
        % but if variables empty : skip
        if isempty(angle)
            disp('to delete : ')
            disp([date fish])
            continue
        end
        
        % --- 2nd criterium : framerate ---
        % .............................................................
        framerate = framerateInfo(:,3);
        lowframerate = find(D.initialization.framerateIni(:,3) < 20);
        if ~isempty(lowframerate)
            del_fHz = lowframerate;
        else
            del_fHz = [];
        end
        
        % --- then check if the experiment hasn't yet been checked
        %     manually ---
        % .............................................................

        checkeddates = find(strcmp(date, EnucSpontUsed(:,1)));
        if ~isempty(checkeddates)
            checkedfish = find(strcmp(fish, EnucSpontUsed(checkeddates,2)));
            cellrow = checkeddates(checkedfish);
            if ~isempty(checkedfish)
                d = EnucSpontUsed(cellrow,3);
                if size(d,1)>1
                    del = cell2mat(d(1));
                else
                    del = cell2mat(d);
                end
            end
            
        else % :: go manual ::
            
            [angleCum] = angle_cum(angle, framerateInfo);
            
            fig = 0;
            angleCumCorr = find_aberrant_bouts(angleCum, coordinates, fig);
            
            fig = 1;
            [boutIndices, ~, ~] = BoutSpot(coordinates, angleCumCorr, framerateInfo, fig);
            
            del = input('sequences to delete manually : ');
            
            if length(del) == size(angleCumCorr,1)
                disp('whole experiment discarded : ')
                disp([fish date])
                continue
                
            end
        end
        
        % --- Finally Delete, if needed, all bad sequences ---
        % .............................................................
        del = del(:);
        del = sort(unique(del));
        if ~isempty(del)
            if size(del,2)> size(del,1)
                del = del';
            end
        end
        DEL = [del_length ; del_fHz ; del];
        DEL = unique(DEL);
        Dini = D.initialization;
        if ~isempty(DEL)
            Dini.angleIni(DEL,:)=[];
            Dini.coordinatesIni(:,:,DEL)=[];
            Dini.framerateIni(DEL,:) = [];
        end
        
        save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
                exptype '/' date '/' fish 'ini.mat'], 'Dini');
        
        % check EXPERIMENT
        %__________________________________________________________________
        
        EnucSpontUsedNew{count, 1} = date;
        EnucSpontUsedNew{count, 2} = fish;
        EnucSpontUsedNew{count, 3} = del;
        
        count = count +1;
        
    end
end

%%
cellcleaning = find(sum(cellfun(@isempty, EnucSpontUsed),2)==3);
EnucSpontUsedNew(cellcleaning,:) = [];
save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/EnucExperimentsUsed.mat'], 'EnucExperimentsUsed'); 
