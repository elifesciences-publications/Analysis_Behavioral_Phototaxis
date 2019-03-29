function[h] = autocorr_function_plots(dXissbiais, FishID, varargin)

% --- mean autocorrelation of dX ---
[ACinorm, sem] = xcorrMatrixRows (dXissbiais);

% fish mean autocorrelation
[ACinormpf, sempf, ~] = xcorrMatrixRowsPerFish (dXissbiais, FishID);

%% plotting
[colour] = colour_palette(0, 1);

bouts = [0:12];

% ***
figure;
ylim([-0.05 1])
ybreaks = [0.2 0.9];

%h = breakerrorbar([bouts' bouts'] ,[ACinorm(bouts+1)' ACinormpf(bouts+1)'], [sem(bouts+1)' sempf(bouts+1)'],...
%    ybreaks(1), ybreaks(2), 'Line');
h = breakerrorbar([bouts'] ,[ACinormpf(bouts+1)'], [sempf(bouts+1)'],...
    ybreaks(1), ybreaks(2), 'Line');
%h=errorbar([bouts'] ,[ACinormpf(bouts+1)'], [sempf(bouts+1)']);
h(1).LineWidth = 2;
h(1).Color = [0.2 0.2 0.2];
% h(2).LineWidth = 2;
% h(2).Color = colour(1,:);

legend('ACF_{fish}');%, '<ACF>_{fish}')
xlim([bouts(1) bouts(end)])
xlabel('bout #')
grid off
title('Autocorrelation function of  \deltaX')
ax = gca;
ax.FontSize = 14;

if nargin > 2
    hold on
    t1 = varargin{1};
    Cor =  varargin{2};
    p1 = breakplot_onexisting(t1(bouts+1),Cor(bouts+1), ybreaks(1), ybreaks(2));
    %p1=plot(t1(bouts+1),Cor(bouts+1));
    p1.DisplayName = 'Analytical solution';
    p1.Color = colour(4,:);
    if nargin > 4
        t2 = varargin{3};
        ac =  varargin{4};
        p2 = breakplot_onexisting(t2(bouts+1), ac(bouts+1), ybreaks(1), ybreaks(2));
        p2.DisplayName = 'Simulation';
        p2.Color = colour_sim(3,:);
    end
end
