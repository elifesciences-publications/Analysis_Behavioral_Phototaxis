function[Es] = spont_save_pooled_exps(Xi, TimeBouti, xCoordi, yCoordi, FishNi)

path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
name = 'spontaneous_swim.mat';

% save structure from spont for 3 experiment types
dxCoord = diff(xCoordi, 1, 2);
dyCoord = diff(yCoordi, 1, 2);
Dist = sqrt(dxCoord.^2 + dyCoord.^2);

%advancement
trajOrientation = wrapToPi(atan(dyCoord./dxCoord));
trajOrientation(dxCoord < 0) = trajOrientation(dxCoord < 0) - pi;
dAlpha = wrapToPi(Xi(:, 1:end-1)) + trajOrientation;
Ri = Dist.*cos(dAlpha);
Ti = Dist.*sin(dAlpha);

Es.Angle = Xi;
Es.FishN = FishNi;
Es.R = Ri;
Es.T = Ti;
Es.TimeBout = TimeBouti;
Es.Coordinates.x = xCoordi;
Es.Coordinates.y = yCoordi;

save([path name], 'Es')

end
