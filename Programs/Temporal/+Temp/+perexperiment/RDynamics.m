
% time dynamics of experiments

%%
colour = Temp.temp_colours;

% mean amplitude
%***
f = figure
plot(nanmean(s6_dx.^2,1), 'Linewidth', 2)
hold on
plot(nanmean(e6_dx.^2,1), 'Linewidth', 2)
plot(nanmean(e3_dx.^2,1), 'Linewidth', 2)
title('mean dx^2 at each bout')
xlim([0 50])
xlabel('bout #')
legend('e60', 'e30', 's60') 

% R and R prj
Rs6 = circ_r(s6_x, [], [], 1);
Re6 = circ_r(e6_x, [], [], 1);
Re3 = circ_r(e3_x, [], [], 1);

ms6 = circ_mean(s6_x, [], 1);
me6 = circ_mean(e6_x, [], 1);
me3 = circ_mean(e3_x, [], 1);

Rps6 = smooth(Rs6.*cos(ms6-pi),3);
Rpe6 = smooth(Re6.*cos(me6-pi),3);
Rpe3 = smooth(Re3.*cos(me3-pi),3);

%***
f = figure;
plot(Rps6, 'Linewidth', 2, 'DisplayName', 'sin60', 'Color', colour(1,:))
hold on
plot(Rpe6, 'Linewidth', 2, 'DisplayName', 'exp60', 'Color', colour(2,:))
plot(Rpe3, 'Linewidth', 2, 'DisplayName', 'exp30', 'Color', colour(3,:))
title('mean resultant vector length (moment) R')
legend
xlim([0 100])
ax = gca;
ax.FontSize = 14;

figure
plot(Rpe6, 'Linewidth', 2)
hold on
plot(Rpe3, 'Linewidth', 2)
plot(Rps6, 'Linewidth', 2)
xlim([0 100])
grid on
ax = gca;
ax.FontSize = 14;

plot(sum(isfinite(e6_x), 1))
hold on
plot(sum(isfinite(e3_x), 1))
plot(sum(isfinite(s6_x), 1))
title('number of bouts')
xlim([0 100])

%***
% figure
% x = e6_x;
% bin = 3;
% for i = 1 : size(x,2) - bin
%     bx = x(:, i : i+bin);
%     polarhistogram(bx(:), 30)
%     text(0, 10, num2str(i))
%     drawnow
%     pause(0.7)
% end

clear me6 me3 ms6