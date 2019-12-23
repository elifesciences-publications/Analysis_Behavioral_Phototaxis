
%% Visualizations for individual fish
uniquefish = unique(FishNi)';
if logical(visu)
    for i = uniquefish(2:end)
        fish = find(FishNi == i);
        %         h = figure;
        %         subplot(3,2,1)
        %         polarhistogram(x(:), 30);
        %         title(i)
        %
        %         subplot(3,2,3)
        %         cmean = circ_mean(x(:,1:50));
        %         r = circ_r(x(:,1:50));
        %         plot(r.*cos(cmean+pi))
        %
        %         subplot(3,2,5)
        %         plot(ones(1, length(fish)), rad2deg(x(:, 1)), '*');
        %         ylim([0 360])
        
        
        % --- look at individual trajectories ---
        visualizeTrajectoryAndAdv(fish, Xi, xCoordi, yCoordi, Disti, Ri, dXi, trajOrientationi)
        hold off
        %         Vart1 = DIlr(fish, 1:end-1);
        %         Vart2 = dX(fish, 1:end);
        %
        %         % --- bins with equal number of elements
        %         [binvals, elts_per_bin, v2bin] = BinsWithEqualNbofElements(Vart1, Vart2, 5, 10);
        %         mV2 = nanmean(v2bin,2);
        %
        %         %***
        %         subplot(3, 2, 6)
        %         errorbar(binvals, mV2, std(v2bin,1,2)/sqrt(elts_per_bin), 'Linewidth', 2)
    end
end
