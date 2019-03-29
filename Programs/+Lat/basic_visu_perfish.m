function[] = basic_visu_perfish(XLat, Xfilt, DIlr, FishN, visux, traj)


%% Visualize trajectories and distributions / fish
% --- individual plots of distributions ---

if logical(visux)
    for i = unique(FishN)'
        fish = find(FishN == i);
        if sum(fish>300) > 0
            break
        end
        x=XLat(fish,:);
        h = figure;
        subplot(3,2,1)
        polarhistogram(x(:), 30);
        title(i)
        
        subplot(3,2,3)
        cmean = circ_mean(x(:,1:50));
        r = circ_r(x(:,1:50));
        plot(r.*cos(cmean+pi))
        
        subplot(3,2,[2])
        histogram(rad2deg(x-Xfilt(fish,:)))

        subplot(3,2,4)
        plot(x');
        
        subplot(3,2,5)
        plot(ones(1, length(fish)), rad2deg(x(:, 1)), '*');
        ylim([0 360])
        
        Vart1 = DIlr(fish, 1:end-1);
        Vart2 = dX(fish, 1:end);
        
        % --- bins with equal number of elements
        [binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, 5, 10);
        mV2 = nanmean(v2bin,2);
        
        %***
        subplot(3, 2, 6)
        errorbar(binvals, mV2, std(v2bin,1,2)/sqrt(elts_per_bin), 'Linewidth', 2)
    end
end

if logical(traj)
    % --- look at individual trajectories ---
    visualizeTrajectoryAndAdv(fish, Xlab, xCoord, yCoord, Dist, Advancement, dXf, trajOrientation, 0)
end

% xcorr = XLat;
%  for i = unique(FishN)'
%          fish = find(FishN == i);
%          x=XLat(fish,:);
% %         cmean = circ_mean(x(:,1:50));
% %         r = circ_r(x(:,1:50));
% %         figure
% %         subplot(1,2,1)
% %         plot(r.*cos(cmean+pi))
% %         subplot(1,2,2)
% %         plot((r.*cos(cmean+pi))/sqrt(length(fish)))
%     fisherr = (x-Xfilt(fish,:));
%     fishmeasureerror(i) = nanmean(fisherr(:));
%     xcorr(fish,:) = XLat(fish,:)-fishmeasureerror(i);
%  end

%% --- distributions per fish --- 
xfn = [];
meanX_perfish = [];
Rperfish = [];
figure;
title('distributions per fish')
for i = sort(unique(FishN),'descend')'
	fish = find(FishN == i);
    xfish = Xfilt(fish,3:20);
    if sum(isnan(xfish(:))) == numel(xfish)
        continue
    end
    polarhistogram(xfish(:),30);
    pause(0.5)
    meanX_perfish(i) = circ_mean(xfish(:));
    Rperfish(i) = circ_r(xfish(:));
    xfn = [xfn; xfish(:)];
end
Rproj_perfish = Rperfish.*cos(meanX_perfish+pi);

%% --- Bias (<X> = f(Il-Ir) ) for individual fish

dX = diff(XLat,1,2);

subplotcount = 1;
%***
figure;

meanXperfish = NaN(length(unique(FishN)),1);
Rcircperfish = NaN(length(unique(FishN)),1);

for i = unique(FishN)'
    fish = find(FishN == i);
    Vart1 = DIlr(fish, 1:end-1);
    Vart2 = dX(fish, 1:end);
    
    % --- bins with equal number of elements
    [binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, 3, 12);
    mV2 = nanmean(v2bin,2);
    
    %***
    subplot(3, 2, subplotcount)
    errorbar(binvals, mV2, std(v2bin,1,2)/sqrt(elts_per_bin), 'Linewidth', 1)
    hold on
    x = Xfilt(fish,:);
    
    %polarhistogram(x(:), 20)
    title([num2str(i) ', seqs ' num2str(length(fish))])
    if mod(subplotcount, 6) == 0
        figure;
        subplotcount = 0;
    end
    subplotcount = subplotcount + 1;
    xfish = XLat(fish,:);
    xfish = xfish(:);
    meanXperfish(i) = circ_mean(xfish);
    Rcircperfish(i) = circ_r(xfish);
end
Rprojperfish = Rcircperfish.*cos(meanXperfish+pi/2);