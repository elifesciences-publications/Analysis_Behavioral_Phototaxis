%% Check distances made at each bout

%% Distribution

dist = Disti(:)/11.5;
dist(isnan(dist)) = [];

%--- pdf distibution ---
N = numel(dist);
binwidth = 0.2;
bins = min(dist):binwidth:max(dist);
[yhist, edges]  = histcounts(dist,bins);
xhist=edges(1:end-1)+binwidth/2;
yhist=yhist/(sum(yhist)*binwidth);

% --- fit ---
% lognormal
distfittype = fittype('1./(x*sigma*sqrt(2*pi)) .* exp(-(log(x)-mu).^2/(2*sigma^2))',...
            'coefficients', {'mu', 'sigma'});
myfit = fit(xhist(2:end)',yhist(2:end)',distfittype, 'startpoint', [1 1]);

%*** 
figure
histogram(dist, 500, 'Normalization', 'pdf', 'FaceAlpha', 0.2, 'DisplayName', 'pdf, small bins')
hold on
plot(xhist, yhist, 'LineWidth', 2, 'DisplayName', 'pdf larger bins')
plot(xhist, myfit(xhist), 'LineWidth', 2, 'DisplayName', 'fit')
legend
xlabel('mm')
ylabel('pdf')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
%% mean/median distance in time

dtimemean = nanmean(Disti,2);
dtimemedian = nanmedian(Disti,2);
plot(dtimemean)
hold on
plot(dtimemedian)

%% autocorrelation in bout distance
dt=20;
var1 = Disti(:,1:end-dt);
var2 = Disti(:,1+dt:end);

% --- randomize for test ---
% var1=var1(:); var1(isnan(var2))=[];
% var2=var2(:); var2(isnan(var2))=[];
% var2 = var2(randperm(numel(var2)))';
% ---

[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, 20, 30);
mv2 = median(v2bin,2);
stdv2 = std(v2bin,1,2);


%*** 
%figure;
errorbar(binvals, mv2, stdv2/sqrt(elts_per_bin))

%% per fish
%***
figure
hold on
for i = unique(FishNi)'
    fish = find(FishNi == i);
    d=Disti(fish,:);
    d=d;%-nanmean(d(:));
    
%     N = numel(d(:))-sum(isnan(d(:)));
%     binwidth = 3*iqr(d(:))/(N^(1/3));
%     bins = min(d(:)):binwidth:max(d(:));
%     [yhist, edges]  = histcounts(d(:),bins)
%     xhist=edges(1:end-1)+binwidth/2;
%     yhist=yhist/(sum(yhist)*binwidth);
%     plot(xhist, yhist)
    
    var1 = d(:,1:end-dt);
    var2 = d(:,1+dt:end);
    
    var1=var1(:); var1(isnan(var2))=[];
    var2=var2(:); var2(isnan(var2))=[];
    var2 = var2(randperm(numel(var2)))';
    
    if size(var1,1) >1
    [binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, 3, 7);
    mv2 = mean(v2bin,2);
    stdv2 = std(v2bin,1,2);

    errorbar(binvals, mv2, stdv2/sqrt(elts_per_bin), 'LineWidth', 1.5, 'Marker', 'o')
    end
end

%% ACF

[ACinorm, sem] = xcorrMatrixRows (Disti);

plot(ACinorm)

%% distance vs dX
var1 = dXi;
var2 = Disti/11.5;

<<<<<<< HEAD
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, 12, 18);
=======
[binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(var1, var2, 20, 24);
>>>>>>> master
mv2 = median(v2bin,2);
stdv2 = std(v2bin,1,2);

%*** 
%figure;
errorbar(binvals, mv2, stdv2/sqrt(elts_per_bin),...
    'LineWidth', 1.5, 'DisplayName', 'median')

mud = NaN(size(v2bin,1),3);
sigmad = NaN(size(v2bin,1),3);
for i = 1 : size(v2bin,1)
    dist = v2bin(i,:);
    binwidth = 0.3;
    bins = min(dist):binwidth:max(dist);
    [yhist, edges]  = histcounts(dist,bins);
    xhist=edges(1:end-1)+binwidth/2;
    yhist=yhist/(sum(yhist)*binwidth);
    
    distfittype = fittype('1./(x*sigma*sqrt(2*pi)) .* exp(-(log(x)-mu).^2/(2*sigma^2))',...
        'coefficients', {'mu', 'sigma'});
    myfit = fit(xhist(2:end)',yhist(2:end)',distfittype, 'startpoint', [1 1]);
    
    figure;
    plot(myfit, xhist, yhist)
    
    cf = confint(myfit)
    mud(i,1) = myfit.mu;
    mud(i,2:3) = cf(:,1)';
    sigmad(i,1) = myfit.sigma;
    sigmad(i,2:3) = cf(:,2)';
end

errorbar(binvals, mud(:,1), mud(:,2)-mud(:,1), mud(:,1)-mud(:,3),...
    'LineWidth', 1.5, 'DisplayName', '\mu')
hold on
errorbar(binvals, sigmad(:,1), sigmad(:,2)-sigmad(:,1), sigmad(:,1)-sigmad(:,3),...
    'LineWidth', 1.5, 'DisplayName', '\sigma')
legend
xlabel('\delta\theta_n')
ylabel('d_n')
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
%%
histogram(Disti(abs(dXi)<0.2))
hold on
histogram(Disti(abs(dXi)>0.2))

