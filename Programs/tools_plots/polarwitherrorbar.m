function [] = polarwitherrorbar(angle,avg,error,fillcolor)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first two input variables ('angle' and 'avg') are same as the input 
% variables for a standard polar plot. The last input variable is the error
% value. Note that the length of the error-bar is twice the error value we
% feed to this function. 
% In order to make sure that the scale of the plot is big enough to
% accommodate all the error bars, i used a 'fake' polar plot and made it
% invisible. It is just a cheap trick. 
% The 'if loop' is for making sure that we dont have negative values  when
% an error value is substrated from its corresponding average value. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_data = length(angle);

fake = polar(angle, 0.12*ones(1,n_data)); %max(avg+error)
set(fake,'Visible','off', 'HandleVisibility', 'off'); 
hold on; 

for ni = 2 : n_data
    if (avg(ni)-error(ni)) < 0
        p = polar(angle(ni)*ones(1,3),[0, avg(ni), avg(ni)+error(ni)]); 
        % change !!!
    else
        p = polar([angle(ni-1)*ones(1,3) angle(ni)*ones(1,3) angle(ni-1)],...
            [avg(ni-1)-error(ni-1), avg(ni-1), avg(ni-1)+error(ni-1) avg(ni)+error(ni), avg(ni),...
            avg(ni)-error(ni), avg(ni-1)-error(ni-1)]); 
        p.LineStyle = 'none';
        p.HandleVisibility = 'off';
    end
    a = fill(get(p,'XData'), get(p,'YData'), fillcolor);
    a.EdgeColor = [1 1 1];
    a.EdgeAlpha = 0;
    a.DisplayName = '<\theta>_f_i_s_h +/- \sigma_f_i_s_h /\surd n';
    if ni ~= n_data
        a.HandleVisibility = 'off';
    end
    set(a,'facealpha',.3)
end

p = polar(angle,avg, 'o');
p.Color = fillcolor;
p.MarkerSize = 3;
p.HandleVisibility = 'off';
legend

hold off
