% diffusion lateralized experiment

Contrasts = DIlr(:,1:23)/max(abs(DIlr(:)));
A = Contrasts*Acoeff;
[MSD] = msdX0shuffled ( XLat(:,1:23) ,A);

% ***
%figure
plot(MSD);
xlim([0 50])
%%
% il faut soustraire le biais à chaque pas de temps

[msd0] = msdX0(XLat, 0);
[MSD] = msdX0shuffled( XLat , 0);

interval = [6:20];
linfit = polyfit(interval, MSD(interval), 1);

%***
figure
plot(interval, MSD(interval))
hold on
plot(interval, linfit(1)*interval+linfit(2))

Dcoeff = linfit(1);
