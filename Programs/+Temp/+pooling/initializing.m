function [L,l] = initializing(Dates, Fish, exptype, intmax, timemin)
% Initializing loop for temporal phototaxis data
% in order to get the size of matrix for pooling data

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
        load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype '/' intmax '/' date '/' fish '/data' fish '.mat']);
        
        % delete the shortest sequences and load the variables of interest
        [~,angle_source,~,framerate,coordinates,~] = remove_sequence(D,timemin); 
        
        % correct the angle
        angle_source_correct = find_aberrant_bouts(angle_source, coordinates, 0);
        
        % get the number of bouts
        [~, ~, nbouts] = BoutSpot(coordinates, angle_source, framerate, 1, 'ini');

        if L < nbouts
            L = nbouts;
        end
        l = l+size(angle_source_correct,1);
    end
end