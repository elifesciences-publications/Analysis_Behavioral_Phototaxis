function[] = visualizeTrajectoryAndAdv(fish, Xlab, xCoord, yCoord, Dist, Adv, dX, trajOrientation, splot)

dxCoord = diff(xCoord, 1, 2);
dyCoord = diff(yCoord, 1, 2);

for s = fish'
    fishorient = (Xlab(s,:));
    dx = dxCoord(s,:);
    dy = dyCoord(s,:);
    d = Dist(s,:);
    adv = Adv(s,:);
    dthe = dX(s,:);
    n = find(adv<0 & abs(dthe)<pi/2);
        
    unitvecty1 = 20*sin(fishorient-pi);
    unitvectx1 = 20*cos(fishorient);
    %
    % trajorient = atan(dy./dx);
    % trajorient(dx < 0) = trajorient(dx < 0) - pi;
    
    trajorient = trajOrientation(s,:);
    unitvecty = d.*sin(trajorient);
    unitvectx = d.*cos(trajorient);
    
    %***
%     if splot
%         subplot(3,2,[2 4])
%     else
%         figure;    
%         subplot(3,1,1:2)
%     end
    plot(xCoord(s,:), yCoord(s,:), '--*', 'Color', [0.8 0.8 0.8])
    hold on
    plot(xCoord(s,1), yCoord(s,1), 'sq')
    plot(xCoord(s,n), yCoord(s,n), 'o')
    hold on
    set(gca,'DataAspectRatio',[1,1,1])
    for i = 1 : length(unitvectx)- sum(isnan(unitvectx))
        uvx = linspace(xCoord(s,i),xCoord(s,i)+unitvectx(i), 10);
        uvy = linspace(yCoord(s,i),yCoord(s,i)+unitvecty(i), 10);
        uvx1 = linspace(xCoord(s,i),xCoord(s,i)+unitvectx1(i), 10);
        uvy1 = linspace(yCoord(s,i),yCoord(s,i)+unitvecty1(i), 10);
        plot(uvx, uvy,'k', uvx1, uvy1, 'r')
    end
    xlim([0 1000])
    ylim([0 1000])
%     if splot
%         subplot(3,2,6)
%     else
%         subplot(3,1,3)
%     end
    %plot(Adv(s,:))
    %title(['seq ' num2str(s) ', ' num2str(cumulativeAdvancement(s))])
    
    waitforbuttonpress
end