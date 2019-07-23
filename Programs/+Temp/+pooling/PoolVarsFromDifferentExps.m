function[E] = PoolVarsFromDifferentExps()

% save structure for 3 experiment types

%% select experiments

Char = {'exp60', 'exp30', 'sin60'};
path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';

for e = 1 : length(Char)
    char = Char{e};
    
    [Dates, Fish, exptype, lum_th, intmax] = Temp.perexperiment.select_exp_type(char);
    
    timemin = 3;
    minBoutNumber = 3;
    
    % initial LOOP
    [L,l] = Temp.pooling.initializing(Dates, Fish, exptype, intmax, timemin);
    
    % pooling LOOP
    [Xrad, Xlab, Lu, TimeBout, xCoord, yCoord, FishN] = Temp.pooling.loop(Dates, Fish, l, L, exptype, intmax, timemin, minBoutNumber);
    
    
    if ( [L, l] - size(Xrad) ) < 0
        warning('problem with array size initial estimation')
    end
    
    dxCoord = diff(xCoord, 1, 2);
    dyCoord = diff(yCoord, 1, 2);
    Dist = sqrt(dxCoord.^2+ dyCoord.^2);

    %advancement
    trajOrientation = wrapToPi(atan(dyCoord./dxCoord));
    trajOrientation(dxCoord < 0) = trajOrientation(dxCoord < 0) - pi;
    dAlpha = wrapToPi(Xlab(:, 1:end-1)) + trajOrientation;
    R = Dist.*cos(dAlpha);
    T = Dist.*sin(dAlpha);
    
    E(e).ExpType = char;
    E(e).AngleSource = Xrad;
    E(e).FishN = FishN;
    E(e).Lum = Lu;
    E(e).Dist = Dist;
    E(e).R = R;
    E(e).TimeBout = TimeBout;
    E(e).T = T;
    E(e).a =  xCoord;
    E(e).b = yCoord;
    E(e).alpha = Xlab;
end

save([path 'temporal_exps.mat'], 'E')

end
