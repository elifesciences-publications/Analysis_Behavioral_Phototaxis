function [Bouts, Sequences] = spont_initial_loop(Dates, Fish, exptype, intmax, timemin)

% shut off polyfit warning
id = 'MATLAB:polyfit:PolyNotUnique';
warning('off',id)

% ini
Bouts=0;
Sequences=0;

for p = 1: length(Dates) % find longest bout sequence
    for q = 1: length(Fish(p,:))
        if iscellstr(Fish(p,q)) == 0
            continue
        end
        % load the needed run in loop
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
        
        % delete the shortest sequences and load the variables of interest
        [angleLab_ini, coordinates_ini, framerate_ini] = remove_sequence_ini(Dini, timemin); 
        
        if isempty(angleLab_ini)
            continue
        end
        
        % transform the angle in lab referential to a cumulated angle
        [angleLabCum_ini] = angle_cum(angleLab_ini, framerate_ini);
        
        % correct the angle
        angleLabCumCorr_ini = find_aberrant_bouts(angleLabCum_ini, coordinates_ini, 0);
        
        % get the number of bouts
        [~, ~, nbouts] = BoutSpot(coordinates_ini, angleLabCumCorr_ini, framerate_ini, 0, 'ini');

        if Bouts < nbouts
            Bouts = nbouts;
        end
        Sequences = Sequences + size(angleLabCumCorr_ini,1);
    end
end