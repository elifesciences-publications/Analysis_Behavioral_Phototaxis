%% plot mean with vector lenght

stat = source_stat_exp30;
vec = stat(:,4);
moy = stat(:,2);

x = vec.*cos(moy*pi/180);
y = vec.*sin(moy*pi/180);

scatter(x,y)
hold on
ang = linspace(0,2*pi);
xcircle = cos(ang);
ycircle = sin(ang);
scatter(xcircle,ycircle,'k')

