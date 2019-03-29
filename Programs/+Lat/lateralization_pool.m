function[El, DatesUsed, FishUsed] = lateralization_pool(path)

%% ALL OK data
Lat.DatesFish

 %% ini loop
[L,l] = Lat.lat_initial_loop(Dates, Fish, exptype, timemin);

%% pooling loop
[XLat, Xlab, Xfilt, TimeBout, xCoord, yCoord, FishN, DatesUsed, FishUsed, CorrespondingFishN] = ...
    Lat.lat_pooling_loop(Dates, Fish, l, L, exptype, timemin, minBoutNumber);

%% save
[El] = Lat.lateralization_save_pooled_data(Xlab, XLat, Xfilt, FishN, TimeBout, xCoord, yCoord,DatesUsed, FishUsed,CorrespondingFishN, path);

end