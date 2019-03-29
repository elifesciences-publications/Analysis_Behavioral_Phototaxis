% variance analysis
% colors = [ 145/255 171/255 60/255;
%     229/255 98/255 38/255;
%     89/255 10/255 49/255];

[colour] = Temp.temp_colours();

colors =  [colour(2,:); colour(3,:); colour(4,:)];


%%
% load variables from E
Temp.perexperiment.loadPooledExpsCell
et{1} = 'sin60';
et{2} = 'exp60';
et{3} = 'exp30';

pxmm = 11.5;

%%
% path for saving figures
path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/PooledData/Figures/';
mkdir(path)

%% stats
for type = 1 : size(et,2)
    
    %choose experiment type
    [data] = Temp.perexperiment.choose_experiment(et{type});
    size(data)
end
%% plain distributions

for type = 1 : size(et,2)
    
    %choose experiment type
   [lum, dlum, x, dx, r, tb, ibi, t, f] = Temp.perexperiment.choose_experiment(et{type});
    color = colors(type,:);
    
    boutsperseq = size(x,2)-sum(isnan(x),2);
    medbout = round(median(boutsperseq));
    
    x = x(:,2:medbout);
    x = x-pi;
    
    Nbins = 18;
    angle = ((1:Nbins) - 0.5) * 2*pi/Nbins;
    xpdf = hist(mod(x(:),2*pi),Nbins)/sum(hist(mod(x(:),2*pi),Nbins));
    
    circr = circ_r(wrapToPi(x(:)));
    theta_mean = circ_mean(wrapToPi(x(:)));
    Rproj_mean = circ_r(x(:)).*cos(theta_mean);
    
    % text on plot
    txtAngle = [pi/2, pi/2];
    txtRng = [0.06, 0.04];
    names = {['<\theta> = ' num2str(round(wrapToPi(theta_mean)*100)/100) ' rad'],...
    ['<R> = ' num2str(round(Rproj_mean*100)/100)]};
    
    %***
    fig1 = figure;
    
    %subplot(2,1,1);
    polarplot([angle, angle(1)],smooth([xpdf, xpdf(1)]),...
        'Linewidth', 2, 'Color', color);
    hold on
    polarplot([0 theta_mean], [0 circr/10], 'Linewidth', 2, 'Color', [0.7 0 0])
    text(txtRng.*cos(txtAngle),txtRng.*sin(txtAngle), names,...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize', 14, 'FontName', 'Times New Roman')
    %polarhistogram(wrapToPi(x(:)),0 : pi/12 :  2*pi,...
    %    'Normalization', 'probability', 'FaceColor', color,'FaceAlpha',.3);
    set(gca,'ThetaZeroLocation','bottom',...
        'ThetaDir','clockwise');
    thetaticks( 0 :90: 360 )
    thetaticklabels({'0' '\pi/2' '\pi' '-\pi/2'})
    title([ et{type} ' : pdf(\theta) bouts 2 to ' num2str(medbout) ', \theta \in [-\pi \pi]'])
    ax = gca;
    ax.RAxisLocation = 45;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    ax.FontWeight = 'normal';
    
    fig3 = polar_distribution(x, Nbins, f, 4);
end

%% plain distributions : ALL ON ONE FIGURE

%***
figure
for type = 1 : size(et,2)
    
    %choose experiment type
   [~, ~, x, ~, ~, ~, ~, ~, f] = Temp.perexperiment.choose_experiment(et{type});
    color = colors(type,:);
    
    boutsperseq = size(x,2)-sum(isnan(x),2);
    medbout = round(median(boutsperseq));
    
    x = x(:,2:medbout);
    x = x-pi;
    
    Nbins = 18;

    polar_distribution_hold_on(x, Nbins, f, colour(type,:));
end


%% R per fish

%***
figure
hold on
for type = 1 : size(et,2)
    
    %choose experiment type
    [lum, dlum, x, dx, r, tb, ibi, t, f] = Temp.perexperiment.choose_experiment(et{type});
    color = colors(type,:);
    meanXperfish = NaN(length(unique(f)),1);
    Rcircperfish = NaN(length(unique(f)),1);
    nbouts = NaN(length(unique(f)),1);
    stdRonseqperfish = NaN(length(unique(f)),1);
    meanRonseqperfish = NaN(length(unique(f)),1);
    seqperfish = NaN(length(unique(f)),1);
    
    medianboutnumber = median(size(x,2)-sum(isnan(x),2));
    
    for i = unique(f)'
        fish = find(f == i);
        xfish = x(fish,2:medianboutnumber)+pi/2;
        meanRonseqperfish(i) = nanmean(circ_r(xfish'));
        stdRonseqperfish(i) = nanstd(circ_r(xfish'));
        seqperfish(i) = length(xfish(:,1)) - sum(isnan(xfish(:,1)));
        xfish = xfish(:);
        xfish(isnan(xfish))=[];
        meanXperfish(i) = circ_mean(xfish);
        Rcircperfish(i) = circ_r(xfish);
        nbouts(i) = numel(xfish);
    end
    Rprojperfish = Rcircperfish.*cos(meanXperfish+pi/2);
    meanRprojonseqperfish = meanRonseqperfish.*cos(meanXperfish+pi/2);
    
    %***
    subplot(2,1,1)
    polarplot(meanXperfish+pi/2, Rcircperfish, 'ko', 'MarkerFaceColor', colour(type,:),...
        'DisplayName', ['Mean on all bouts/individual in ' et{type}])
    hold on
    set(gca,'ThetaZeroLocation','bottom',...
        'ThetaDir','clockwise')
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    ax.ThetaAxisUnits = 'radians';
    ax.RAxisLocation = pi;
    ax.RMinorGrid = 'on';
    rticks(0:0.5:1)
    thetaticks((0 : pi/4 : 2*pi))
    thetaticklabels({'0', '', '\pi/2', '', '\pi', '', '-\pi/2'})
    legend
    
    subplot(2,1,2)
    plot(wrapToPi(meanXperfish+pi/2), Rprojperfish,...
        'ko', 'MarkerFaceColor', colour(type,:), 'HandleVisibility', 'off')
    hold on
    errorbar(wrapToPi(meanXperfish+pi/2),meanRprojonseqperfish, stdRonseqperfish./sqrt(seqperfish),...
        'Color', colour(type,:), 'Marker', '.', 'MarkerEdgeColor', colour(type,:), 'LineStyle', 'none',...
        'DisplayName', ['mean +/- sem in ' et{type}])
    plot([-pi pi], [0 0], '--k', 'HandleVisibility', 'off');
    plot([0 0], [-1 1], '--k', 'HandleVisibility', 'off');
    
    xlim([-pi pi])
    xticks([-pi,-pi/2, 0,pi/2, pi])
    xticklabels({'-\pi','-\pi/2', '0','\pi/2', '\pi'})
    xlabel('<\theta>_f_i_s_h')
    ylabel('R_f_i_s_h')
    ax=gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 14;
    legend
end

%% ........................................................................
% BINS WITH UNEVEN NUMBER OF EVENTS

% choose variables to plot and dt
v1 = 3;
dt = 1;

%***
fig = figure;
r = rand;
scatcolor = [0 r r];
for type = 1 : size(et,2)
    %choose experiment type
    [lum, dlum, x, dx, r, tb, ibi, t, f] = Temp.perexperiment.choose_experiment(et{type});
    color = colors(type,:);
    
    if v1 == 1           % memory of previous dX
        Vart1 = dlum(:, 1:end-dt);
        Vart2 = dx(:, dt+1:end).*(dx(:, 1:end-dt)./abs(dx(:, 1:end-dt)));
        xl = 'dLu';
        yl = 'dX(n)x dX(n-1)/|dX(n-1)|';
    elseif v1 == 2       % r vs dlum
        Vart1 = dlum(:, 1:end-dt)./((lum(:,1:end-dt-1)+lum(:,1+dt:end-dt))/2);
        Vart2 = r(:, dt+1:end)/pxmm;
        xl = 'dLu/Lu';
        yl = 'R (mm)';
    elseif v1 == 2.5       % r vs dlum
        Vart1 = dlum(:, 1:end-dt);
        Vart2 = r(:, dt+1:end)/pxmm;
        xl = 'dLu';
        yl = 'R (mm)';
    elseif v1 == 3       % ibi vs dlum
        Vart1 = dlum(:, 1:end-dt)./((lum(:,1:end-dt-1)+lum(:,1+dt:end-dt))/2);
        Vart2 = ibi(:, dt+1:end);
        xl = 'dLu';
        yl = 'inter-bout interval';
    elseif v1 == 4       % dx^2 vs dlum
        Vart1 = dlum(:, 1:end-1)./((lum(:,1:end-2)+lum(:,2:end-1))/2);
        Vart2 = dx(:, 2:end).^2;
        xl = 'dLu/Lu';
        yl = 'dX^2';
    elseif v1 == 5       % dx^2 vs dlum
        Vart1 = dlum(:, 1:end-1);%./((lum(:,1:end-2)+lum(:,2:end-1))/2);
        Vart2 = t(:, dt+1:end)/pxmm;
        xl = 'dLu/Lu';
        yl = 'dX^2';
    end
    
    b1 = prctile(Vart1(:), 99.9);
    b2 = prctile(Vart1(:), 0.1);
    bins = linspace(b2, b1, 30);
    w = mean(diff(bins));
    
    meanV2 = NaN(length(bins), 1);
    varV2 = NaN(length(bins), 1);
    semV2 = NaN(length(bins), 1);
    for i = 2 : length(bins)
        selectedV2 = Vart2(Vart1 <= bins(i) & Vart1 >= bins(i-1));
        meanV2(i-1) = nanmean(selectedV2);
        varV2(i-1) = nanvar(selectedV2);
        n = length(selectedV2);
        semV2(i-1) = nanstd(selectedV2)/sqrt(n);
    end
    
    %***
    if v1 == 2
        yyaxis left
    end
    scat = scatter(Vart1(:), Vart2(:), 7, 'DisplayName', et{type});
    scat.MarkerEdgeAlpha = 0.05;
    scat.MarkerEdgeColor = color;
    hold on
    e = errorbar(bins + w/2, meanV2, semV2, 'DisplayName', et{type});
    e.Color = color;
    e.LineWidth = 1.5;
    ylim([min(Vart2(:)) prctile(Vart2(:), 99)])
    xlim([bins(1) bins(end)])
    xlabel(xl)
    ylabel(yl)
    
    if v1 == 2
        yyaxis right
        plot(bins + w/2, varV2, 'sq', 'MarkerSize', 7, 'MarkerFaceColor', 'r')
        ylabel('R^2')
    end
    hold on
end
grid on
title(['dt = ' num2str(dt)])
legend('show')
ax = gca;
ax.FontSize = 14;

%% ........................................................................
% BINS WITH SAME NUMBER OF ELEMENTS

[colour] = Temp.temp_colours('dark');

% choose variables
v1 = 3;
v2 = 2;
dt = 1;

%***
fig = figure;

minbin = inf;
maxbin = -inf;
miny = inf;
maxy = -inf;
for type = 1 : size(et,2)
    
    %choose experiment type
    [lum, dlum, x, dx, r, tb, ibi, t, f] = Temp.perexperiment.choose_experiment(et{type});
    
    fishdxmean=[];
    for i = unique(f)'
        dxf = dx(i==f,:);
        fishdxmean = [fishdxmean ; repmat(nanmean(dxf(:)), [size(dxf,1) 1])];
    end
        
    
    color = colour(type,:);
    if v1 == 1
        Vart1 = lum(:, 1:end-(1+dt));
        xl = 'I_n (mW.cm^-^2)';
        xticks(0:0.05:0.2)
    elseif v1 == 2
        Vart1 = dlum(:, 1:end-dt);
        xl = '\deltaI_n_-_1';
    elseif v1 == 3
        Vart1 = dlum(:, 1:end-dt)./((lum(:,1:end-dt-1)+lum(:,1+dt:end-dt))/2);
        xl = '\deltaI/I_{n-1} = 2(I_{n-1} - I_{n-2})/(I_{n-1} + I_{n-2})';
    elseif v1 == 4
        Vart1 = dx(:, 1:end-dt);
        xl = '\delta\theta_n_-_1)';    
    elseif v1 == 5
        Vart1 = wrapToPi(x(:, 1:end-dt));
        xl = '\theta_n'; 
        xticks([-pi -pi/2 0 pi/2 pi])
        xticklabels({'-\pi' '-\pi/2' '0' '\pi/2' '\pi'})
    end
    
    if v2 == 1 %dx
        Vart2 = dx(:, 1:end) - fishdxmean;
        yl = '<\delta\theta_n> (rad)';
    elseif v2 == 2 %dx^2
        Vart2 = dx(:, 1+dt:end).^2;
        Vart2(Vart2>100) = nan;
        yl = '<\delta\theta_n^2> (rad)';
    elseif v2 == 3
        Vart2 = dx(:, dt+1:end).*(dx(:, 1:end-dt)./abs(dx(:, 1:end-dt)));
        yl = '\delta\theta_n.\delta\theta_n_-_1/|\delta\theta_n_-_1|';
    elseif v2 == 4       % R
        Vart2 = r(:, dt+1:end)/pxmm;
        yl = 'R(mm)';
    elseif v2 == 5       % ibi
        Vart2 = ibi(:, dt+1:end);
        yl = 'IBI(s)';
    elseif v2 == 6       % transverse
        Vart2 = t(:, dt+1:end)/pxmm;
        yl = 'transverse displacement';
    end
    
    [binvals, elts_per_bin, v2binMatrix] = BinsWithEqualNbofElements(Vart1, Vart2, 12, 17);        
    
    %***
    %scat = scatter(Vart1(:), Vart2(:), 3, 'MarkerEdgeColor', color, 'DisplayName', et{type});
    %scat.MarkerEdgeAlpha = 0.1;
    hold on
    errorbar(binvals, mean(v2binMatrix,2), std(v2binMatrix,1,2)/sqrt(elts_per_bin),...
        'Linewidth', 2, 'Color', color, 'CapSize', 3, 'DisplayName', et{type})
    
    minbin = min(minbin, binvals(1));
    maxbin = max(maxbin, binvals(end));
    miny = min(miny, prctile(Vart2(:), 5));
    maxy = max(maxy, prctile(Vart2(:), 95));
end
xlabel(xl)
ylabel(yl)
ylim([miny maxy])
xlim([minbin*0.8 maxbin*1.2])

grid on
title(['dt = ' num2str(dt)])
legend('show')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%% dX = f(X)

%***
fig = figure;
colors = [ 145/255 171/255 60/255; 
    229/255 98/255 38/255; 
    89/255 10/255 49/255];
minbin = inf;
maxbin = -inf;
miny = inf;
maxy = -inf;
for type = 1 : size(et,2)
    
    %choose experiment type
    [lum, dlum, x, dx, r, tb, ibi, t, f] = Temp.perexperiment.choose_experiment(et{type});
    color = colors(type,:);
    
    xw = wrapToPi(x(:, 1:end-2));
    Vart1_a = xw(xw>0 & xw<pi);
    Vart1_b = xw(xw<0 & xw>-pi);
    xl = 'X_n_-_1';
    xticks([-pi -pi/2 0 pi/2 pi])
    xticklabels({'-\pi' '-\pi/2' '0' '\pi/2' '\pi'})
    
    Vart2_a = dx(:,2:end);% - nanmean(dx,2);
    Vart2_a = Vart2_a(xw>0 & xw<pi);
    Vart2_a(abs(Vart2_a)>pi) = NaN;
    Vart2_b = dx(:, 2:end);% - nanmean(dx,2);
    Vart2_b = Vart2_b(xw<0 & xw>-pi);
    yl = '<dX_n> (rad)';
    
    [binvals_a, elts_per_bin_a, v2bin_a] = BinsWithEqualNbofElements(Vart1_a, Vart2_a, 15, 17);
    [binvals_b, elts_per_bin_b, v2bin_b] = BinsWithEqualNbofElements(Vart1_b, Vart2_b, 15, 17);
    
    %***
    hold on
    scat = scatter(Vart1_a(:), Vart2_a(:), 3, 'MarkerEdgeColor', color, 'DisplayName', et{type});
    scat.MarkerEdgeAlpha = 0.1;
    scat = scatter(Vart1_b(:), Vart2_b(:), 3, 'MarkerEdgeColor', color, 'DisplayName', et{type});
    scat.MarkerEdgeAlpha = 0.1;
    errorbar(binvals_a, mean(v2bin_a,2), std(v2bin_a,1,2)/sqrt(elts_per_bin_a),...
        '-', 'Linewidth', 2, 'Color', color, 'CapSize', 3, 'DisplayName', et{type})
    errorbar(binvals_b, mean(v2bin_b,2), std(v2bin_b,1,2)/sqrt(elts_per_bin_b),...
        '--', 'Linewidth', 2, 'Color', color, 'CapSize', 3, 'DisplayName', et{type})
    minbin = min(minbin, binvals_b(1));
    maxbin = max(maxbin, binvals_a(end));
    miny = min(miny, prctile(Vart2_a(:), 1));
    maxy = max(maxy, prctile(Vart2_a(:), 95));
end
xlabel(xl)
ylabel(yl)
ylim([miny maxy])
xlim([minbin*1.2 maxbin*1.2])

grid on
title(['dt = ' num2str(dt)])
legend('show')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

%% check after a certain number of bouts (shift)

%choose experiment type
%--------------------------------------------------------------------------
type = 3;

[lum, dlum, x, dx, r, tb, ibi, t] = Temp.perexperiment.choose_experiment(et{type});

% chose the 2 variables to plot one against eachother
% -------------------------------------------------------------------------
m1 = 0;
m2 = 30;
figure
for i = m1 : 5 : m2
    shift = i;
    
    Vart1 = dlum(:, 1 :end - (1 + shift))./( (lum(:,2 :end - (1 + shift)) + lum(:,1 :end- (2+shift)))./2 ) ;
    Vart2 = dx(:, 2+shift:end).^2;
    [binvals, elts_per_bin, v2binMatrix] = BinsWithEqualNbofElements(Vart1, Vart2, 20, 50);

    %***
    %f = figure;
    %plot(binvals, mean(v2sqbinMatrix,2), '*')
    hold on
    prettycolor = [i/(m2+20) i/(m2+20) i/(m2+20)];
    errorbar(binvals, mean(v2binMatrix,2), std(v2binMatrix,1,2)/sqrt(elts_per_bin), 'Color', prettycolor)
end

xlabel('dlum(n-1)./( (lum(n-2) +lum(n-1))./2 )')
ylabel('mean dx^2')
title([et{type} ' - shift : ' num2str(m1) ' to ' num2str(m2)])
%% fit dX distribution with 2 gaussians
% -------------------------------------------------------------------------

%choose experiment type
%--------------------------------------------------------------------------

ets = sort(et);

% ***
fig = figure;
for type = 1 : 3
    [lum, dlum, x, dx, r, tb, ibi, t] = Temp.perexperiment.choose_experiment(ets{type});
    
    % variables
    Vart1 = dlum(:, 1 :end-1)./( (lum(:,2:end-1) + lum(:,1:end-2))./2 ) ;
    Vart2 = dx(:, 2:end);
    [binvals, elts_per_bin, v2binMatrix] = BinsWithEqualNbofElements(Vart1, Vart2, 30, 50);
    
    %
    f1 = NaN(size(v2binMatrix,1), 3);
    f2  = NaN(size(v2binMatrix,1), 3);
    nf = NaN(size(v2binMatrix,1), 1);
    nt = NaN(size(v2binMatrix,1), 1);
    for i = 1 : size(v2binMatrix,1)
        
        [dx_histogram, x_histogram] = hist(v2binMatrix(i,:), 2*round(sqrt(elts_per_bin)));
        f = fit(x_histogram.',dx_histogram.','gauss2');
        forward = @(x) f.a1*exp(-((x-f.b1)/f.c1).^2);
        side = @(x) f.a2*exp(-((x-f.b2)/f.c2).^2);
        %     figure;
        %     plot(x_histogram, dx_histogram)
        %     hold on
        %     plot(x_histogram, f.a1*exp(-((x_histogram-f.b1)/f.c1).^2) + f.a2*exp(-((x_histogram-f.b2)/f.c2).^2))
        f1(i, 1) = f.a1;
        f1(i, 2) = f.b1;
        f1(i, 3) = f.c1;
        f2(i, 1) = f.a2;
        f2(i, 2) = f.b2;
        f2(i, 3) = f.c2;
        nf(i) = sum(forward(x_histogram));
        nt(i) = sum(side(x_histogram));
    end
    
    % plot all fits
    % ***
    % figure
    % for i = 1 : size(v2binMatrix,1)
    %     plot(binvals, f1(i,1)*exp(-((binvals-f1(i,2))/f1(i,3)).^2) + f2(i,1)*exp(-((binvals-f2(i,2))/f2(i,3)).^2))
    %     hold on
    % end
    
    c = [rand rand rand];
%     subplot(1,4,1)
%     plot(binvals, f1(:,1), 'o', 'Color', c)
%     hold on
%     plot(binvals, f2(:,1), 'sq', 'Color', c, 'MarkerFaceColor', c)
%     title('amplitude')
%     ax = gca;
%     ax.FontSize = 14;
%     
%     subplot(1,4,2)
%     plot(binvals, f1(:,2), 'o', 'Color', c)
%     hold on
%     plot(binvals, f2(:,2), 'sq', 'Color', c, 'MarkerFaceColor', c)
%     title('mean')
%     ax = gca;
%     ax.FontSize = 14;
%     
%     subplot(1,4,3)
%     plot(binvals, f1(:,3), 'o', 'Color', c)
%     hold on
%     plot(binvals, f2(:,3), 'sq', 'Color', c, 'MarkerFaceColor', c)
%     title('std')
%     ax = gca;
%     ax.FontSize = 14;
%     
%     subplot(1,4,4)
    plot(binvals, nt./nf, 'o', 'Color', c)
    
    hold on
    
end
legend(sort([et et]))

%% keep only std

for type = 1 : 3
    [lum, dlum, x, dx, r, tb, ibi, t] = Temp.perexperiment.choose_experiment(et{type});
    
    % variables
    Vart1 = dlum(:, 1 :end-1)./( (lum(:,2:end-1) + lum(:,1:end-2))./2 ) ;
    Vart2 = dx(:, 2:end);
    [binvals, elts_per_bin, v2binMatrix] = BinsWithEqualNbofElements(Vart1, Vart2, 28, 28);
    
    %
    for i = 1 : size(v2binMatrix,1)
        [dx_histogram, x_histogram] = hist(v2binMatrix(i,:), 2*round(sqrt(elts_per_bin)));
        f = fit(x_histogram.',dx_histogram.','gauss2');
        ft(type, i) = f.c2;
    end
    b(type, :) = binvals;
end
plot(b',ft', 'o')

%***
figure;
mb = mean(b);
mft = mean(ft);
plot(mb, mft, 'o')
p = polyfit(mb(1:13), mft(1:13), 1)
hold on
plot(mb(1:13), mb(1:13)*p(1)+p(2))

%% evolution in time
% load experiment
Temp.perexperiment.loadPooledExpsCell

% time dyn script
Temp.perexperiment.RDynamics

