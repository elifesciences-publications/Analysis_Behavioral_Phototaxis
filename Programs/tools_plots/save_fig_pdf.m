function save_fig_pdf(fig, folder, name)

root_path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Figures/';

if ~exist([root_path folder], 'dir')
    mkdir([root_path folder])
end

set( fig, 'Units','Inches');
pos = get( fig, 'Position');
set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(fig, [root_path folder filesep name], '-dpdf', '-r0')
disp(['figure' name ' saved in :'])
disp([root_path folder])

end