function[Ee, DatesUsed, FishUsed] = enuc_pool(path)

%% ALL OK data
[Dates, Fish, exptype] = Enuc.DatesFish();
timemin = 5;
minBoutNumber = 3;

 %% ini loop
[L,l] = Enuc.initial_loop(Dates, Fish, exptype, timemin);

%% pooling loop
[Xsource, Xlab, Xfilt, TimeBout, xCoord, yCoord, FishN, DatesUsed, FishUsed, Luminosity] = ...
    Enuc.pooling_loop(Dates, Fish, l, L, exptype, timemin, minBoutNumber);

%% save
[Ee] = Enuc.save_pooled_data(Xlab, Xsource, Xfilt, FishN, TimeBout, xCoord, yCoord, Luminosity, path);

end