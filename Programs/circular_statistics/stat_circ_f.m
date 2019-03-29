function [stat] = stat_circ_f(data,fig)

%% Comments
% Do circular statistics for one fish ie use AngleLum structure
%Inputs:
% -load Perfish
%Outputs:
% -Circ_stat_lab: matrix. Each line represents the statistics of one fish.
% 1st column: nombre of bout
% 2nd column: circular mean
% 3rd column: standart deviation
% 4th column: mean vector length
% 5th column: p-value of the rayleight test
% 6th column: p-value for a directionnality test toward the light (180°)
% -Figure

%% make statistics
data = mod(data,360);
stat = nan(size(data,1),6);
for fish = 1:size(data,1)
    f = find(isnan(data(fish,:)),1);
    if isempty(f)==1
        f=size(data,2)+1;
    end
    a = data(fish,1:f-1)'*pi/180;
    nb_bout = size(a,1);
    moy = mod(circ_mean(a)*180/pi,360);
    std = circ_std(a);
    vector = circ_r(a);
    p_ray = circ_rtest(a);
    p180 = circ_vtest(a,pi);
    
    if fig == 1
        figure
        subplot(1,2,1)
        rose(a)
        subplot(1,2,2)
        histogram(a*180/pi,20)
        xlim([0 360])
    end
    
    stat(fish,:) = [nb_bout,moy,std,vector,p_ray,p180];
end