<<<<<<< HEAD
function[] = visualizeTrajectoryAndAdv(fish, X, xCoord, yCoord, Dist, Adv, dX, trajOrientation)

% pix to mm
xCoord = xCoord/11.5;
yCoord = yCoord/11.5;
Dist = Dist/11.5;
Adv = Adv/11.5;

% differentiation
dxCoord = diff(xCoord, 1, 2);
dyCoord = diff(yCoord, 1, 2);

start_color = 'k';
traject_color = [0.5 0 0.5];

for s = fish'
    fishorient = X(s,:);
    dx = dxCoord(s,:);
    dy = dyCoord(s,:);
    d = Dist(s,:); % for display
    adv = Adv(s,:);
    dthe = dX(s,:);
    norm_turn_magn = (dthe+pi)/(2*pi);
    norm_turn_magn = norm_turn_magn - min(norm_turn_magn);
    norm_turn_magn = norm_turn_magn/max(abs(norm_turn_magn));
    n = find(adv<0 & abs(dthe)<pi/2);
        
    unitvecty1 = 2.5*sin(fishorient-pi);
    unitvectx1 = 2.5*cos(fishorient);
=======
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
>>>>>>> master
    %
    % trajorient = atan(dy./dx);
    % trajorient(dx < 0) = trajorient(dx < 0) - pi;
    
    trajorient = trajOrientation(s,:);
    unitvecty = d.*sin(trajorient);
    unitvectx = d.*cos(trajorient);
    
<<<<<<< HEAD
    plot(xCoord(s,:), yCoord(s,:), '-o',...
        'Color', traject_color, 'Markersize', 2, 'MarkerFaceColor', 'k')
    hold on
    plot(xCoord(s,1), yCoord(s,1), 'sq',...
        'Markersize', 7, 'MarkerFaceColor', start_color, 'MarkerEdgeColor', start_color)
=======
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
>>>>>>> master
    plot(xCoord(s,n), yCoord(s,n), 'o')
    hold on
    set(gca,'DataAspectRatio',[1,1,1])
    for i = 1 : length(unitvectx)- sum(isnan(unitvectx))
        uvx = linspace(xCoord(s,i),xCoord(s,i)+unitvectx(i), 10);
        uvy = linspace(yCoord(s,i),yCoord(s,i)+unitvecty(i), 10);
        uvx1 = linspace(xCoord(s,i),xCoord(s,i)+unitvectx1(i), 10);
        uvy1 = linspace(yCoord(s,i),yCoord(s,i)+unitvecty1(i), 10);
        plot(uvx, uvy,'k', uvx1, uvy1, 'r')
<<<<<<< HEAD
        arrow_color = [norm_turn_magn(i) 0 1-norm_turn_magn(i)];
        if isnan(arrow_color)
            return
        end
        drawArrow([uvx1(1) uvx1(end)], [uvy1(1) uvy1(end)], arrow_color);
    end
    
    % plot the ROI contour
    theta = 0:0.05:2*pi+0.05;
    xc = 40;
    yc = 40;
    rc = 40;
    hold on;
    plot(xc+rc*cos(theta),yc+rc*sin(theta), 'k', 'Linewidth', 2);
    
    xlim([0 80])
    ylim([0 80])
    
    xlabel('mm')
    ylabel('mm')
    
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'TimesNewRoman';
    ax.TickLength = [0 0];

    disp(s)
    waitforbuttonpress
    hold off
=======
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
>>>>>>> master
end