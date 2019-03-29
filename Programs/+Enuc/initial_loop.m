function [L,l] = initial_loop(Dates, Fish, exptype, timemin)

% shut off polyfit warning
id = 'MATLAB:polyfit:PolyNotUnique';
warning('off',id)

L=0;
l=0;
 
for p = 1: length(Dates) % find longest bout sequence
    for q = 1: length(Fish(p,:))
        if iscellstr(Fish(p,q)) == 0
            continue
        end
        % load the needed run in loop
        fish = char(Fish(p,q));        
        date = char(Dates(p));
        
        datapath = ['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype filesep date filesep fish '.mat'];
        if ~exist(datapath, 'file')
            continue
        end
        load(datapath, 'E');
        
        disp([ date fish])
        
        % delete the shortest sequences and load the variables of interest
        [~,angle_source,~,framerate,coordinates,~] = remove_sequence(E,timemin); 
        
        if isempty(angle_source)
            continue
        end
        % correct the angle
        angle_source_correct = find_aberrant_bouts(angle_source, coordinates, 0);
        
        % get the number of bouts
        [~, ~, nbouts] = BoutSpot(coordinates, angle_source, framerate, 0, 'ini');

        if L < nbouts
            L = nbouts;
        end
        l = l+size(angle_source_correct,1);
    end
end