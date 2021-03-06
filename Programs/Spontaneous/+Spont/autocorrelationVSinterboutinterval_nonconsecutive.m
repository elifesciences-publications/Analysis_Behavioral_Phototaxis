function[f1, f2, binvals1, mv1, stdv1] = autocorrelationVSinterboutinterval_nonconsecutive(dX, interboutint, turn_thresh, dt, label, varargin)

% autocorrelation of successive bouts as a function of interbout interval
% and extraction of pflip

%% Autocorrelation not only consecutive bouts : up to dt

dX(abs(dX)<turn_thresh) = 0;
% dX(dX<0) = -1;
% dX(dX>0) = 1;

var1 = NaN(dt*size(dX,1), size(dX,2)-1);
var2 = NaN(dt*size(dX,1), size(dX,2)-1);

timebout = [zeros(size(interboutint,1),1), cumsum(interboutint,2)];

for i = 1 : dt
    dXt = dX(:,1:end-i);
    dXtpdt = dX(:,1+i:end);
    pdXndXnp1 = dXt.*dXtpdt;
    pabsdXndXnp1 = abs(dXt).*abs(dXtpdt);
    
    % plug into variables
    var1((i-1)*size(dX,1)+1:(i-1)*size(dX,1)+size(dX,1), 1:end-i+1) = timebout(:,1+i:end-1) - timebout(:,1:end-i-1);
    var2((i-1)*size(dX,1)+1:(i-1)*size(dX,1)+size(dX,1), 1:end-i+1) = pdXndXnp1./pabsdXndXnp1;
end

bins = 18;

[binvals1, elts_per_bin1, v2bin1] = BinsWithEqualNbofElements(var1, var2, bins, bins+3);
mv1 = mean(v2bin1,2);
stdv1 = std(v2bin1,1,2);

%% pflip 
% pturn_estimate = abs(dXt);
% flips = 0.5*(1 - pdXndXnp1./pabsdXndXnp1./pturn_estimate);
% [binvals2, elts_per_bin2, v2bin2] = BinsWithEqualNbofElements(var1, flips, bins, bins+3);
% mv2 = mean(v2bin2,2);
% stdv2 = std(v2bin2,1,2);

%% plots

%***
if nargin > 5
    f1 = varargin{1};
    set(0, 'currentfigure', f1)
    hold on
else
    f1 = figure;
end
hold on
shadedErrorBar(binvals1, mv1, stdv1/sqrt(elts_per_bin1),...
    'lineprops',{'Linewidth', 2, 'DisplayName', [label 'turns only |d\theta|>' num2str(turn_thresh) ', bins : ' num2str(bins)]})
xlabel('inter-bout interval (sec)')
ylabel('\delta\theta_n\delta\theta_{n+1}/|\delta\theta_n||\delta\theta_{n+1}|')


legend
ax=gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 14;

%*** 
if nargin > 5
    f2 = varargin{2};
    set(0, 'currentfigure', f2)
    hold on
else
    f2 = figure;
end

% shadedErrorBar(binvals2, mv2, stdv2/sqrt(elts_per_bin2),...
%      'lineprops',{'LineWidth', 2, 'DisplayName', [label 'turns only |d\theta|>0.22rad']})
% legend