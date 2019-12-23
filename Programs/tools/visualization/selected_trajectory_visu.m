load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
    'lateralisation/18-10-16/fish2.mat'], 'E')

theta_t = E.experiment.angleCum;

seq_of_interest = 5;

traj_of_interest = theta_t(seq_of_interest,:);
traj_of_interest(isnan(traj_of_interest))=[];

fHz = E.experiment.framerate(seq_of_interest,3);
time = [1:length(traj_of_interest)]/fHz;

plot(time, deg2rad(traj_of_interest))
yticks([-pi/2:pi/2 : 2*pi])
yticklabels({'-\pi/2','0','\pi/2', '\pi','3\pi/2', '2\pi'})
xlabel('time (s)')
ylabel('\theta')

ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 18;

ax.Children.LineWidth = 2;