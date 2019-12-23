function[h] = autocorr_function_plots(dXissbiais, FishID, plot_type, varargin)

% --- mean autocorrelation of dX ---
[ACinorm, sem] = xcorrMatrixRows (dXissbiais);

% fish mean autocorrelation
[ACinormpf, sempf, ~] = xcorrMatrixRowsPerFish (dXissbiais, FishID);

%% plotting
[colour] = colour_palette(0, 1);

bouts = 0:12;

% ***
figure;
ylim([-0.05 1])
ybreaks = [0.2 0.9];

switch plot_type
    case 'bp'
    h = breakerrorbar(bouts' ,ACinormpf(bouts+1)', sempf(bouts+1)',...
    ybreaks(1), ybreaks(2), 'Line');
    case 'np'
    h = errorbar(bouts, ACinormpf(bouts+1), sempf(bouts+1));
    hold on
end
h(1).LineWidth = 2;
h(1).Color = [0.2 0.2 0.2];

legend('ACF_{fish}');
xlim([bouts(1) bouts(end)])
xlabel('bout #')
grid off
title('Autocorrelation function of  \deltaX')
ax = gca;
ax.FontSize = 14;

if nargin > 3
    hold on
    t1 = varargin{1};
    Cor =  varargin{2};
    switch plot_type
        case 'bp'
        p1 = breakplot_onexisting(t1(bouts+1), Cor(bouts+1), ybreaks(1), ybreaks(2));
        case 'np'
    	p1 = plot(t1(bouts+1),Cor(bouts+1));
    end
    p1.DisplayName = 'Analytical solution';
    p1.Color = colour(4,:);
    p1.LineWidth = 1.5;
    if nargin > 5
        t2 = varargin{3};
        ac =  varargin{4};
        switch plot_type
            case 'bp'
                p2 = breakplot_onexisting(t2(bouts+1), ac(bouts+1), ybreaks(1), ybreaks(2));
            case 'np'
                p2 = plot(t2(bouts+1), ac(bouts+1));
        end
        p2.DisplayName = 'Simulation';
        p2.Color = colour(3,:);
        p2.LineWidth = 1.5;
    end
end
