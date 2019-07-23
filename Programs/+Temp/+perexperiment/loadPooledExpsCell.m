
%%

path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/';
file = 'temporal_exps';
load([path file], 'E')

%%
for i = 1 : size(E, 2)
    if E(i).ExpType == 'exp60'
        e6_x = E(i).AngleSource;
        e6_lum = E(i).Lum;
        e6_f = E(i).FishN;
        e6_r = E(i).R;
        e6_t = E(i).T;
        e6_tb = E(i).TimeBout;
        e6_a = E(i).a;
        e6_b = E(i).b;
        e6_alpha = E(i).alpha;
    elseif E(i).ExpType == 'exp30'        
        e3_x = E(i).AngleSource;
        e3_lum = E(i).Lum;
        e3_f = E(i).FishN;
        e3_r = E(i).R;
        e3_t = E(i).T;
        e3_tb = E(i).TimeBout;
        e3_a = E(i).a;
        e3_b = E(i).b;
        e3_alpha = E(i).alpha;
    elseif E(i).ExpType == 'sin60'
        s6_x = E(i).AngleSource;
        s6_lum = E(i).Lum;
        s6_f = E(i).FishN;
        s6_r = E(i).R;
        s6_t = E(i).T;
        s6_tb = E(i).TimeBout;
        s6_a = E(i).a;
        s6_b = E(i).b;
        s6_alpha = E(i).alpha;
    end
end

%%
% luminosity profiles
e3_profile = luminosity_exponentielle(0.3);
e3_p_angle = deg2rad(e3_profile(:,1));
e3_p_int = e3_profile(:,2);

e6_profile = luminosity_exponentielle(0.6);
e6_p_angle = deg2rad(e6_profile(:,1));
e6_p_int = e6_profile(:,2);

s6_profile = luminosity_sinus(0.6);
s6_p_angle = deg2rad(s6_profile(:,1));
s6_p_int = s6_profile(:,2);

%%
% other useful variables
% number of fish per experiment type
e3_N = length(unique(e3_f));
e6_N = length(unique(e6_f));
s6_N = length(unique(s6_f));

% dlum
e3_dlum = diff( e3_lum, 1, 2 );
e6_dlum = diff( e6_lum, 1, 2 );
s6_dlum = diff( s6_lum, 1, 2 );
% dx
e3_dx =  diff( e3_x, 1, 2 );
e6_dx =  diff( e6_x, 1, 2 );
s6_dx =  diff( s6_x, 1, 2 );

%%
clear E