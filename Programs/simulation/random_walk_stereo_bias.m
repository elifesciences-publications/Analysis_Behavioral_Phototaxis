function[theta_complete, lL_R] = random_walk_stereo_bias(Nexp, Ntimes, psw_turn, p_turn, wturn, wstraight, mturn, varargin)

% INPUT
% Nexp number of experiments
% Ntimes number of time steps per experiment

% --- initialization ---
rng('shuffle')
theta_complete = NaN(Nexp, Ntimes); % complete dataset
lL_R = NaN(Nexp, Ntimes);

if nargin > 7
    initial_bias = varargin{8};
    xdata = varargin{9};
    if logical(initial_bias)
        theta_ini = randsample(xdata(:,1), Nexp, true); % bias the initial distribution
    end
else
    theta_ini = rand(1,Nexp)*2*pi; % initial angles uniformly distributed
end

lum_lr_ini = interp1([-pi 0 pi], [1 -1 1], wrapToPi(theta_ini));

lrst = (rand > 0.5)*2-1; % intial turn state random 50/50

% --- walking loop ---
tic
for N = 1:Nexp
    theta = NaN(1,Ntimes);
    lum_lr = NaN(1,Ntimes);
    theta(1) = theta_ini(N);
    lum_lr(1) = lum_lr_ini(N);
    for t = 1:Ntimes-1
        th = theta(t);
        lrst = lrst*((rand > psw_turn)*2-1); % independant from whole field illum
        turn = rand <= p_turn;
        if turn
            dth = lrst * (wturn*randn) - mturn(lum_lr(t));
        else
            dth = lrst * (wstraight*randn);
        end
        theta(t+1) = th + dth;
        lum_lr(t+1) = interp1([0 pi], [-1 1], abs(wrapToPi(theta(t+1))));
    end
    theta_complete(N,:) = theta;
    lL_R(N,:) = lum_lr;
end
