function[theta, luminosity, lambda] = random_walk_temporal_bias(Nexp, Ntimes, lum, theta_ini)

% --- base parameters ---
pflip = 0.19;       % probabibility of switching left vs right states
p_turn = 0.41;      % probability of triggering a turn swim (otherwise go straight)

wturn=0.59;         % standard deviation of turning distribution
wstraight=0.092;    % 0.12

% --- intial state ---
%theta_ini = rand(Nexp,1)*2*pi - pi; % initial angles uniformly distributed
luminosity_prev = 50*ones(Nexp,1); % luminosity during adaptation
luminosity_ini = interp1(lum(:,1), lum(:,2), abs(wrapTo180(rad2deg(theta_ini))));
lrst = (rand(Nexp,1) > 0.5)*2-1; % intial turn state random 50/50
turn_last = (rand(Nexp,1) > 0.5)*2-1;

% --- initialize variables ---
theta = [theta_ini NaN(Nexp, Ntimes-1)]; % complete dataset
luminosity = [luminosity_ini NaN(Nexp, Ntimes-1)];
lambda = NaN(Nexp, Ntimes-1);

tic
for t = 1:Ntimes
    % --- variable from previsou cycle ---
    theta_last = theta(:,t);
   
    lrst = lrst.*((rand(Nexp,1) > pflip)*2-1); % independent from whole field illum
    
    if t == 1
        dll = 8*(luminosity_ini - luminosity_prev)./((luminosity_ini + luminosity_prev)/2); 
    else
        dll = (luminosity(:,t) - luminosity(:,t-1))./((luminosity(:,t) + luminosity(:,t-1))/2) ;
    end
    [p_turn, wturn] = ProbaTurnFLum(dll);
    
    %aplimtude corr
    %p_turn(turn_last == 0) = 0.8*p_turn(turn_last == 0);
    %p_turn(turn_last == 1) = 1.1*p_turn(turn_last == 1);
    
    % --- separate turns and forward swims ---
    % default : forward scoots
    dth = wstraight*randn(Nexp,1);  % delta theta
    l = gamrnd(2.5, 1, Nexp,1);     % distance
    
    % select turns
    
    turn = rand(Nexp,1) < p_turn;
    turn_last = turn; % store
    dth(turn == 1) = lrst(turn == 1) .* abs(wturn(turn==1).*randn(sum(turn==1),1));
    l(turn==1) = gamrnd(3.3, 1, sum(turn),1);
    
    % --- set new values ---
    theta(:,t+1) = theta_last + dth; 
    luminosity(:,t+1) = interp1(lum(:,1), lum(:,2), abs(wrapTo180(rad2deg(theta(:,t+1)))));
    lambda(:,t) = l;
end
toc
