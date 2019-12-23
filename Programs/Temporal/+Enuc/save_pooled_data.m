function[Ee] = save_pooled_data(Xlab, XLat, Xfilt, FishN, TimeBout, xCoord, yCoord, Luminosity, path)

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

Ee.ExpType = exptype;
Ee.AngleSource = XLat;
Ee.AngleLab = Xlab;
Ee.AngleSourceFiltered = Xfilt;
Ee.Luminosity = Luminosity;
Ee.FishN = FishN;
Ee.xCoord = xCoord;
Ee.yCoord = yCoord;
Ee.R = R;
Ee.TimeBout = TimeBout;
Ee.T = T;

save([path 'enucleated_exps.mat'], 'Ee')
