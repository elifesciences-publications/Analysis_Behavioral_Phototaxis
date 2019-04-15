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

%% Amplitude decorrelation

fig = Spont.Amplitude(dXi);

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