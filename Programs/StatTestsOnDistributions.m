
% Statistical tests for all data

%% Load data & simulations
Spont.Load;
x0 = wrapToPi(Xi);

Lat.Load;
xl = wrapToPi(XLat - 3*pi/2); % fix the source on 0

Temp.Load
xt = wrapToPi(X - pi);

clearvars -except x0 xt xl

[~, ~, xt_e6, ~, ~, ~, ~, ~, ~,] = Temp.perexperiment.choose_experiment('exp60');
[~, ~, xt_e3, ~, ~, ~, ~, ~, ~,] = Temp.perexperiment.choose_experiment('exp30');
[~, ~, xt_s6, ~, ~, ~, ~, ~, ~,] = Temp.perexperiment.choose_experiment('sin60');

xt_e6 = wrapToPi(xt_e6 - pi);
xt_e3 = wrapToPi(xt_e3 - pi);
xt_s6 = wrapToPi(xt_s6 - pi);

typei = 'spont';
typel = 'lat';
typet = 'temp';
typet_e6 = 'temp_e6';
typet_e3 = 'temp_e3';
typet_s6 = 'temp_s6';

%% R-test on non-uniformity of data
% If pval -> 0 ==> non-uniform distribution
[p0, z0] = circ_rtest(x0(:));
[pl, zl] = circ_rtest(xl(:));
[pt, zt] = circ_rtest(xt(:));

[pt_e6, zt_e6] = circ_rtest(xt_e6(:));
[pt_e3, zt_e3] = circ_rtest(xt_e3(:));
[pt_s6, zt_s6] = circ_rtest(xt_s6(:));

%% V test for non-uniformity of circular data with a specified mean
% direction dir. ---
% H0: the population is uniformly distributed around the circle
% pval -> 0 ==> rejection of H0
% HA: the population is not distributed uniformly around the circle but has a mean of dir.
%
%   Note: Not rejecting H0 may mean that the population is uniformly
%   distributed around the circle OR that it has a mode but that this mode
%   is not centered at dir.

dir_source = 0;
[pv0, v0] = circ_vtest(x0(:), dir_source);

[pvl, vl] = circ_vtest(xl(:), dir_source);

[pvt, vt] = circ_vtest(xt(:), dir_source);

[pvt_e6, vt_e6] = circ_vtest(xt_e6(:), dir_source);
[pvt_e3, vt_e3] = circ_vtest(xt_e3(:), dir_source);
[pvt_s6, vt_s6] = circ_vtest(xt_s6(:), dir_source);


%% Table output
tests = {'pval'; 'zscore'; 'pval'; 'vscore'};

type = {typei, typel, typet, typet_e6, typet_e3, typet_s6};
Pz = {p0, pl , pt, pt_e6, pt_e3, pt_s6};
Z = {z0 , zl, zt, zt_e6, zt_e3, zt_s6};
Pv = {pv0, pvl, pvt, pvt_e6, pvt_e3, pvt_s6};
V = {v0, vl, vt, vt_e6, vt_e3, vt_s6};

varNames = ['Tests', type];
varTypes = {'string', 'double', 'double', 'double', 'double', 'double', 'double'};
sz = [size(tests,1) size(varNames,2)];
stat_table = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames',varNames);
stat_table(:,1) = tests;
stat_table(1,2:end) = Pz;
stat_table(2,2:end) = Z;
stat_table(3,2:end) = Pv;
stat_table(4,2:end) = V;

disp(stat_table)

%% Perform simulations

% --- general parameters ---
psw_turn = 0.19; 
p_turn = 0.41; 
wturn= 0.6;
wstraight= 0.09;
mturn = @(x) -pi/6*x; 

Nexp = 1000;
Ntimes = 30;

% --- lateralized ---
[theta_complete, lL_R] = random_walk_stereo_bias(Nexp, Ntimes, psw_turn, p_turn, wturn, wstraight, mturn);

xsl = theta_complete - 3*pi/2;
typesl = 'simu_lat';

% --- temporal ---
lum_exp6 = luminosity_exponentielle(0.6);
lum_exp3 = luminosity_exponentielle(0.3);
lum_sin6 = luminosity_sinus(0.6);

[theta_complete, ~] = random_walk_temporal_bias(Nexp, Ntimes, psw_turn, p_turn, wturn, wstraight, lum_exp6);
xst_e6 = theta_complete;
typest_e6 = 'simu_temp_e6';

[theta_complete, ~] = random_walk_temporal_bias(Nexp, Ntimes, psw_turn, p_turn, wturn, wstraight, lum_exp3);
xst_e3 = theta_complete;
typest_e3 = 'simu_temp_e3';

[theta_complete, ~] = random_walk_temporal_bias(Nexp, Ntimes, psw_turn, p_turn, wturn, wstraight, lum_sin6);
xst_s6 = theta_complete;
typest_s6 = 'simu_temp_s6';

%% Test ?? 
% Parametric Watson-Williams multi-sample test for equal means : not
% possible because the resultant vector length has to be R>0.45

[p] = two_sample_distribution_test(xl, xsl)


[pwt_e6, k, K] = circ_kuipertest(test(:), test2(:), 100, 1)
[pwt_e3, k, K] = circ_kuipertest(xt_e3_lin(:), xst_e3_lin(:), 18, 1)
[pwt_s6, ~] = circ_wwtest(xt_s6(:), xst_s6(:))

