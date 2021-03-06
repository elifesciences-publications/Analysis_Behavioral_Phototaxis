%% Analysis of spontaneous swim sequences (constant whole field luminous stimulation)

%% ************************************************************************

% --- Load data ---
Spont.Load

%% Individual fish raw visualization
visu = 1;

Spont.IndividualFishVisu

%% General stats

% --- number of sequences VS number of bouts
NofSeqs_ini = size(Xi,1);
totalNumberOfBouts_ini = length( Xi(:) ) - sum( isnan( Xi(:) ) );
boutsperseq = sort( size(Xi,2) - sum(isnan(Xi),2), 'descend');

% ***
figure
subplot(2,2,[1:2])
plot( boutsperseq, 1 : NofSeqs_ini )
xlabel('number of bouts')
ylabel('number of sequences')
title('Number of bouts per sequence')

subplot(2,2,3)
histogram(TimeBouti)
text(150, 2000, { ['max : ' num2str(nanmax(TimeBouti(:)))]});
xlabel('sec')
title('time of bout')

subplot(2,2,4)
histogram(IBIi)
text(5, 2000, {['mean : ' num2str(nanmean(IBIi(:)))], ['median : ' num2str(nanmedian(IBIi(:)))]});
xlim([0 10])
xlabel('sec')
title('interbout interval')

%% X distributions 
% ***
Nbins= 18;
fig = polar_distribution(Xi, Nbins, FishID, 1);
fig.Name = 'Spontaneous';

%% trajectories

xCoord0 = (xCoordi-xCoordi(:,1))/pxmm;
yCoord0 = (yCoordi-yCoordi(:,1))/pxmm;

randomization = repmat(-pi*rand(size(xCoord0,1),1),1,size(xCoord0,2));
x = xCoord0.*cos(randomization) - yCoord0.*sin(randomization);
y = xCoord0.*sin(randomization) + yCoord0.*cos(randomization);

figure
plot(x',y','.-')
title('no stimulation')
xlabel('mm')
ylabel('mm')

%% --- scatter of wrapped dX & mean/var on bins ---
edit Spont.dXvsX

%% %% dX distribution
[Wturn, Wfor, Pturn] = Spont.dXdistribution(dXi, different_fish)
wturn = Wturn(1);
wfor = Wfor(1);
pturn = Pturn(1);

%% <dX> = A (Bias,Drift) per bout for each sequence / each fish ---
save_fig = 0;
Spont.bias_distribution_individuals...
    (different_fish, dXi, FishID, sequencesperfish, fishdXmean, fishdXstd, save_fig)

%% Amplitude correlation

fig = Spont.Amplitude(dXi);

% autocorrelation of dXsq
dT = 10; % lag
Spont.autocorr_diff_lag_dt(dXissbiais.^2,dT);

[ACsqinormpf, semsqpf, acPF] = xcorrMatrixRowsPerFish ((dXissbiais.^2), FishID);

errorbar([0:30], ACsqinormpf(1:31), semsqpf(1:31))

%% autocorrelation of dXsq per fish 
ind_fish = unique(FishID);
dT = 1;

mV1 = nan (length(ind_fish), 20);
binvals = nan (length(ind_fish), 20);
for i = 1 : length(ind_fish)
    seq_fish = (FishID==ind_fish(i));
    dxseq = dXissbiais(seq_fish,:);
    if numel(dxseq) - sum(isnan(dxseq(:))) < 50
        disp('d�gag�')
        continue
    end
    [a, b] = Spont.autocorr_diff_lag_dt_square(dxseq,dT);
    mV1(i, 1:length(a)) = a;
    binvals(i,1:length(b)) = b;
end

% mean on all fish :
amplitude_corr_meanfish = mV1;
allbinvals = sort(binvals(:));
allbinvals(isnan(allbinvals)) = [];

interp_mV1 = NaN(length(ind_fish), length(allbinvals));
for i = 1 : length(ind_fish)
    a = binvals(i,:);
    a(isnan(a)) = [];
    b = mV1(i,:);
    b(isnan(b)) = [];
    if isempty(a)
        continue
    end
    interp_mV1(i,:) = interp1(a, b, allbinvals);
end

%random permutation
interp_mV1_rand = nan(size(interp_mV1));
for i = 1 : length(ind_fish)
    interp_mV1_rand(i,:) = interp_mV1(i,randperm(length(allbinvals)));
end

%***
%figure
% % randomized
% plot(allbinvals, interp_mV1_rand, '.', 'MarkerSize', 1.5, 'Color', [0.9 0.9 0.9])
% hold on
% % real
% plot(allbinvals, interp_mV1, '.', 'MarkerSize', 1.5)
% randomized mean
shadedErrorBar(allbinvals, nanmean(interp_mV1_rand), nanstd(interp_mV1_rand)/sqrt(length(ind_fish)-3),...
    'lineprops',{'.','Color', [0.5 0.5 0.5]})
%real mean
shadedErrorBar(allbinvals, nanmean(interp_mV1), nanstd(interp_mV1)/sqrt(length(ind_fish)-3),...
    'lineprops',{'-','Color', [0.2 0.4 0.9], 'LineWidth', 1.5})

ylabel('<\delta\theta^2_{n+1}>')
xlabel('<\delta\theta^2_{n}>')
legend('<\delta\theta^2_{randomized}>', '<\delta\theta^2_{n+1}>')

ax = gca;
ax.XScale = 'log';
ax.FontName = 'TimesNewRoman';
ax.FontSize = 14;


%% Autocorrelation analysis
edit Spont.Autocorrelation

edit Spont.AutocorrelationTime

%% Diffusion analysis
edit Spont.Diffusion

%% Look at advancement (distance between 2 bouts)

edit Spont.Distance

%% Inter-bout interval
histogram(IBIi)
ibi = IBIi;
binwidth = 0.1;
bins = min(ibi(:)):binwidth:max(ibi(:));
[yhist, edges]  = histcounts(ibi,bins);
xhist=edges(1:end-1)+binwidth/2;
yhist=yhist/(sum(yhist)*binwidth);

distfittype = fittype('1./(x*sigma*sqrt(2*pi)) .* exp(-(log(x)-mu).^2/(2*sigma^2))',...
    'coefficients', {'mu', 'sigma'});
myfit = fit(xhist(2:end)',yhist(2:end)',distfittype, 'startpoint', [1 1]);

figure;
plot(myfit, xhist, yhist)

% ---