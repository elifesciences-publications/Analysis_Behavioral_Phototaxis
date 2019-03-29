function[El] = lateralization_save_pooled_data...
    (Xlab, XLat, Xfilt, FishN, TimeBout, xCoord, yCoord, DatesUsed, FishUsed, CorrespondingFishN, path)

exptype = 'lateralization';

dxCoord = diff(xCoord, 1, 2);
dyCoord = diff(yCoord, 1, 2);
Dist = sqrt(dxCoord.^2+ dyCoord.^2);

%advancement
trajOrientation = wrapToPi(atan(dyCoord./dxCoord));
trajOrientation(dxCoord < 0) = trajOrientation(dxCoord < 0) - pi;
dAlpha = wrapToPi(Xlab(:, 1:end-1)) + trajOrientation;
R = Dist.*cos(dAlpha);
T = Dist.*sin(dAlpha);

El.ExpType = exptype;
El.AngleSource = XLat;
El.AngleLab = Xlab;
El.AngleSourceFiltered = Xfilt;
El.FishN = FishN;
El.xCoord = xCoord;
El.yCoord = yCoord;
El.R = R;
El.TimeBout = TimeBout;
El.T = T;
Linfo = {DatesUsed; FishUsed; CorrespondingFishN};

save([path 'lateralized_exps.mat'], 'El')
save([path 'lateralized_info.mat'], 'Linfo')
