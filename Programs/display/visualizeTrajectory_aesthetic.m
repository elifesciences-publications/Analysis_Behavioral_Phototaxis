function[] = visualizeTrajectory_aesthetic(fish, Xlab, xCoord, yCoord, Dist, Adv, dX, trajOrientation)

s = fish;

fishorient = (Xlab(s,:));
d = Dist(s,:);
r = Adv(s,:);
dthe = dX(s,:);

unitvecty1 = 2*sin(fishorient-pi);
unitvectx1 = 2*cos(fishorient);
%
% trajorient = atan(dy./dx);
% trajorient(dx < 0) = trajorient(dx < 0) - pi;

trajorient = trajOrientation(s,:);
unitvecty = d.*sin(trajorient);
unitvectx = d.*cos(trajorient);

%***
plot(xCoord(s,:), yCoord(s,:),...
    '--o', 'Color', [0.2 0.2 0.3], 'MarkerSize', 3, 'MarkerFaceColor', [0 0 0])
hold on
plot(xCoord(s,1), yCoord(s,1), 'sq')
hold on
set(gca,'DataAspectRatio',[1,1,1])
for i = 1 : length(unitvectx)- sum(isnan(unitvectx))+1
    % between bouts
    uvx = linspace(xCoord(s,i),xCoord(s,i)+unitvectx(i), 10);
    uvy = linspace(yCoord(s,i),yCoord(s,i)+unitvecty(i), 10);
    % orientation vector
    uvx1 = linspace(xCoord(s,i),xCoord(s,i)+unitvectx1(i), 10);
    uvy1 = linspace(yCoord(s,i),yCoord(s,i)+unitvecty1(i), 10);
    plot(uvx, uvy,'k'... , uvx1, uvy1, 'r')
        )
    % arrows
    drawArrow([uvx1(1) uvx1(end)], [uvy1(1) uvy1(end)])
end

xlim([0 100])
ylim([0 100])
