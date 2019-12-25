% Load simulated data from Seb

load('/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/DataSimulation/withEC/sequence_delta_teta.mat')
load('/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/DataSimulation/withEC/time_bout.mat')
load('/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/DataSimulation/withEC/type_of_bout.mat')
dXsimuEC = Sequence_delta_teta(1:end-1);
ibisimuEC = diff(time_bout);
typeboutsimuEC = type_spikeL_all_bout;
load('/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/DataSimulation/withoutEC/sequence_delta_teta.mat')
load('/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/DataSimulation/withoutEC/time_bout.mat')
load('/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/DataSimulation/withoutEC/type_of_bout.mat')
dXsimunoEC = Sequence_delta_teta(1:end-1);
ibisimunoEC = diff(time_bout);
typeboutsimunoEC = type_spikeL_all_bout;

% --- global turn threshold ---
turn_thresh = 0.22;

% --- with EC ---
dX = dXsimuEC;
interboutint = ibisimuEC/10;

Spont.autocorrelationVSinterboutinterval(dX, interboutint, turn_thresh, 'simu EC', fig1, fig2)


% --- without EC ---
dX = dXsimunoEC;
interboutint = ibisimunoEC/10;

Spont.autocorrelationVSinterboutinterval(dX, interboutint, turn_thresh, 'simu no EC', fig1, fig2)
