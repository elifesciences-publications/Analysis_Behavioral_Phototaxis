function [] = save_fig_and_svg(fig, folder, name)

root_path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/Figures/';

if ~exist([root_path folder], 'dir')
    mkdir([root_path folder])
end
    
saveas(fig,[root_path folder filesep name],'fig')
saveas(fig,[root_path  folder filesep  name],'svg')
disp('fig saved .fig & .svg in :')
disp([root_path folder])

end