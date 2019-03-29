%% Comments
% Do circular statistics for one fish ie use AngleLum structure
%Inputs:
% -load AngleLum
%Outputs:
% -Circ_stat_lab: matrix. Each line represents the statistics of one fish.
% 1st column: 



%% Code
clearvars -except ang_bout_r ang_bout_cum 
close all

ang_lab = one_sequence(ang_bout_r)'*pi/180;
ang_source = one_sequence(ang_bout_cum)'*pi/180;

%% ----- For the laboratory angle -----
m_ang_lab = mod(circ_mean(ang_lab)*180/pi,360);
r_v_mean_lab = circ_r(ang_lab);
p_ray_lab = circ_rtest(ang_lab);

if p_ray_lab > 0.05
    disp('p-value rayleight lab test > 0.05')
else
    disp('p-value raileight lab test < 0.05')
    p180lab = circ_vtest(ang_lab, 180*pi/180);
end

for i=0:360
    plab(i+1) = circ_vtest(ang_lab, i*pi/180);
end
plot(plab)

%% ----- For the source angle -----
m_ang_source = mod(circ_mean(ang_source)*180/pi,360);
r_v_mean_source = circ_r(ang_source);
p_ray_s = circ_rtest(ang_source);

if p_ray_s > 0.05
    disp('p-value rayleight source test > 0.05')
else
    disp('p-value raileight source test < 0.05')
    p180s = circ_vtest(ang_source, 180*pi/180);
end

for i=0:360
    ps(i+1) = circ_vtest(ang_source, i*pi/180);
end
figure
plot(ps)