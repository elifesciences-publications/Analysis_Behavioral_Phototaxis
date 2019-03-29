
nData = size(data,1); %# number of data points
numberOfdeltaT = floor(nData/4); %# for MSD, dt should be up to 1/4 of number of data points

msd = zeros(numberOfDeltaT,3); %# We'll store [mean, std, n]

%# calculate msd for all deltaT's

for dt = 1:numberOfDeltaT
   deltaCoords = data(1+dt:end,2:4) - data(1:end-dt,2:4);
   squaredDisplacement = sum(deltaCoords.^2,2); %# dx^2+dy^2+dz^2

   msd(dt,1) = mean(squaredDisplacement); %# average
   msd(dt,2) = std(squaredDisplacement); %# std
   msd(dt,3) = length(squaredDisplacement); %# n
end