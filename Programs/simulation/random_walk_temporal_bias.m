function[theta, luminosity, lambda] = random_walk_temporal_bias(Nexp, Ntimes, lum)

% --- base parameters ---
psw_turn = 0.119; % probabibility of switching left vs right states
p_turn = 0.41; % probability of triggering a turn swim (otherwise go straight)

wturn=0.59; % 0.74
wstraight=0.092; % 0.12
% ---

theta_ini = rand(Nexp,1)*2*pi; % initial angles uniformly distributed
luminosity_ini = interp1(lum(:,1), lum(:,2), abs(wrapTo180(rad2deg(theta_ini))));
lrst = (rand(Nexp,1) > 0.5)*2-1; % intial turn state random 50/50

theta = [theta_ini NaN(Nexp, Ntimes-1)]; % complete dataset
luminosity = [luminosity_ini NaN(Nexp, Ntimes-1)];
lambda = NaN(Nexp, Ntimes-1);

tic
for t = 1:Ntimes
    theta_last = theta(:,t);
    lrst = lrst.*((rand(Nexp,1) > psw_turn)*2-1); % independent from whole field illum
    if t > 2
        dll = (luminosity(:,t) - luminosity(:,t-1))./((luminosity(:,t) + luminosity(:,t-1))/2) ;
        [p_turn, wturn] = ProbaTurnFLum(dll);
    end
    % --- separate turns and forward swims ---
    turn = rand(Nexp,1) < p_turn;
    % when turns
    dth = lrst .* abs(wturn.*randn(Nexp,1)) .* turn;
    %l = turn .* pick_from_lognorm_distrib(possible_lambda_values, sigma_dturn, muturn, Nexp)';
    l = gamrnd(2.8, 1, Nexp,1);
    l(turn==1) = gamrnd(3.5, 1, sum(turn),1);
    
%     % when scoots
%     if sum(turn==0)>0
%         dth(turn==0) = wstraight*randn(sum(turn==0),1);
%         l(turn==0) =  pick_from_lognorm_distrib(possible_lambda_values, sigma_dfwd, mufwd, sum(turn==0))';
%     end
    
    % --- set new values ---
    theta(:,t+1) = theta_last + dth; 
    luminosity(:,t+1) = interp1(lum(:,1), lum(:,2), abs(wrapTo180(rad2deg(theta(:,t+1)))));
    lambda(:,t) = l;
    
    if mod(t,10) == 0
        disp(t)
    end
end
