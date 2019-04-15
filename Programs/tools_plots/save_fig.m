function [] = save_fig(fig, folder, name, format)

root_path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/AnalysisOutput/Figures/';

if ~exist([root_path folder], 'dir')
    mkdir([root_path folder])
end

set( fig, 'Units','Inches');
pos = get( fig, 'Position');
set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
    
switch format
    case {'fig', 'svg'}
        saveas(fig,[root_path folder filesep name], format)
        
    case {'png', 'pdf'}
        print(fig, [root_path folder filesep name], ['-d' format], '-r0')
        
    case 'all'
        % fig
        saveas(fig,[root_path folder filesep name], 'fig')
        
        % svg
        saveas(fig,[root_path folder filesep name], 'svg')
        
        % pdf
        print(fig, [root_path folder filesep name], '-dpdf', '-r0')
        
        %png
        print(fig, [root_path folder filesep name], '-dpng', '-r0')
        
end

disp(['figure' name '.' format ' saved in :'])
disp([root_path folder])

end