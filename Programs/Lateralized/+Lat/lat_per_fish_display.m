
fish = find(FishID == 5);
xfish = XLat(fish, bornes(1):bornes(2) )+pi/2;

circr = circ_r(wrapToPi(xfish(:)))
theta_mean = circ_mean(wrapToPi(xfish(:)))


%***
polarscatter(xfish(:), ones(size(xfish(:))), 'ok')
hold on
polarplot([0 theta_mean], [0 circr])

ax = gca;
ax.ThetaZeroLocation = 'bottom';
ax.Children(1).LineWidth = 1.5;
legend('data points', ['R = ' num2str(circr) ', <\theta> = ' num2str(theta_mean)])
