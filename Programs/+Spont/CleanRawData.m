% clean D before analysis

%% Spont
timemin = 5;

%% load previously stored parameters
load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/SpontExperimentsUsed.mat'], 'SpontExperimentsUsedNew');
SpontExperimentsUsed = SpontExperimentsUsedNew;
clear SpontExperimentsUsedNew

%%
%SpontExperimentsUsedNew = cell(50,3);
%count = 1;
N=4;

L = 0 ;
l = 0 ;
Larg = NaN(N,1);
long = NaN(N,1);
for EXP = 3:N
    if EXP == 1
        [Dates, Fish, exptype, ~, intmax] = Temp.chooseExpType('exp60');
    elseif EXP == 2
        [Dates, Fish, exptype, ~, intmax] = Temp.chooseExpType('exp30');
    elseif EXP ==3
        [Dates, Fish, exptype, ~, intmax] = Temp.chooseExpType('sin60');
    elseif EXP == 4
        Lat.DatesFish;
        intmax = 0;
    end
    for p = 29: length(Dates)
        for q = 1 : length(Fish(p,:))
            
            if iscellstr(Fish(p,q)) == 0
                continue
            end
            fish = char(Fish(p,q));
            date = char(Dates(p));
            
            if EXP < 4
                load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
                exptype '/' intmax '/' date '/' fish '/data' fish '.mat'], 'D');
            else
                load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
                exptype '/' date '/' fish '/data' fish '.mat'], 'D');
            end
            disp([date ' ' fish])
            
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
%            checkeddates = find(strcmp(date, SpontExperimentsUsed(:,1)));
%             if ~isempty(checkeddates)
%                 checkedfish = find(strcmp(fish, SpontExperimentsUsed(checkeddates,2)));
%                 cellrow = checkeddates(checkedfish);
%                 if ~isempty(checkedfish)
%                     d = SpontExperimentsUsed(cellrow,3);
%                     if size(d,1)>1
%                         del = cell2mat(d(1));
%                     else
%                         del = cell2mat(d);
%                     end
%                 
%                 
%                 else % :: go manual ::
                    
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
        %end
             %end
            
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
            
            SpontExperimentsUsedNew{count, 1} = date;
            SpontExperimentsUsedNew{count, 2} = fish;
            SpontExperimentsUsedNew{count, 3} = del;
            
            count = count +1;
            
            if intmax
                save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
                exptype '/' intmax '/' date '/' fish 'ini.mat'], 'Dini');
            else
                save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
                exptype '/' date '/' fish 'ini.mat'], 'Dini');
            end
            
        end
    end
end
%%
cellcleaning = find(sum(cellfun(@isempty, SpontExperimentsUsedNew),2)==3);
SpontExperimentsUsedNew(cellcleaning,:) = [];
save(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/SpontExperimentsUsed.mat'], 'SpontExperimentsUsedNew'); 
