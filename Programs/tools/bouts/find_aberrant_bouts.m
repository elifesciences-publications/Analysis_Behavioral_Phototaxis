function [angle_source_correct] = find_aberrant_bouts(angle_source, coordinates, fig)

% corrects aberrant bouts / punctual or cumulated head-tail revearsal
% spotted dtheta(n) > maxheadingAngle (virtual, just to check : approx. 180-3*std(dtheta))
%
% - delete peak in last point (if there is) in theta
%
% - simple head/tail inversion peak (about 180? peak) :
%   dtheta(n) = -dtheta(n+1)
%
% - head/tail inversion but for 2 consecutive points
%   dtheta(n) = dtheta(n+2)
%   dtheta(n-1) = dtheta(n+1) (around 0)
%
% - head/tail inversion step (180?)
%   dtheta(n-1) = dtheta(n+1) (without overshoot)
%
% - double cumulated head/tail inversion (360?)
%   dtheta(n) = dtheta(n+1)

%%
coordinates(coordinates==0) = NaN;

x = squeeze(coordinates(:,1,:))';
y = squeeze(coordinates(:,2,:))';

maxheadingAngle = 150;
dtheta = diff(angle_source, 1,2);
dthetastd = nanmean(nanstd(dtheta,1,2));
dxy = diff(x, 1, 2).^2 + diff(y, 1, 2).^2;
dxystd = nanmean(nanstd(dxy,1,2));
[misrow, miscol] = find(abs(dtheta) > maxheadingAngle);

angle_source_correct = angle_source;
sopt = 6;

for b = 1 : length(misrow)
    if miscol(b) <= 1 || miscol(b)
        break
    elseif miscol(b) > 1 && miscol(b) <= sopt
        s = miscol(b) - 1;
    else 
        s = sopt;
    end
    dtn = dtheta(misrow(b), miscol(b));
    dtn1 = dtheta(misrow(b), miscol(b)+1);
    if size(dtheta,2) > miscol(b) +1
        dtn2 = dtheta(misrow(b), miscol(b)+2);
    else
        dtn2 = [];
    end
    dtn_1 = dtheta(misrow(b), miscol(b)-1);
    if find(isnan(dtheta(misrow(b),:)), 1)-1 == miscol(b) 
        % last val
        angle_source_correct(misrow(b), miscol(b)+1) = angle_source_correct(misrow(b), miscol(b));
        dtheta(misrow(b), miscol(b)) = 0;
    elseif abs( dtn + dtn1 ) < 2*dthetastd 
        % simple peak
        angle_source_correct(misrow(b), miscol(b)+1) = mode(angle_source_correct(misrow(b), miscol(b)-s : miscol(b)));
        dtheta(misrow(b), miscol(b)) = 0;
        dtheta(misrow(b), miscol(b)+1) = 0;
    elseif abs( dtn_1 - dtn1 ) < 2*dthetastd 
        % either a double peak or a simple step
        if (abs( dtn + dtn2 ) < 2*dthetastd) && (dxy(misrow(b), miscol(b)) < 3*dxystd)
            %double peak
            angle_source_correct(misrow(b), miscol(b)+1) = mode(angle_source_correct(misrow(b), miscol(b)-s : miscol(b)));
            angle_source_correct(misrow(b), miscol(b)+2) = mode(angle_source_correct(misrow(b), miscol(b)-(s-1) : miscol(b)+1));
            dtheta(misrow(b), miscol(b)) = 0;
            dtheta(misrow(b), miscol(b)+1) = 0;
            dtheta(misrow(b), miscol(b)+2) = 0;
        else % simple step
            angle_source_correct(misrow(b), miscol(b)+1:end) = angle_source_correct(misrow(b), miscol(b)+1:end) - dtheta(misrow(b), miscol(b));
            dtheta(misrow(b), miscol(b)) = 0;
        end
        dtheta(misrow(b), miscol(b)) = 0;
    elseif (abs( dtn - dtn1 ) < 2*dthetastd) && (dxy(misrow(b), miscol(b)) < 6*dxystd)
        % double step
        angle_source_correct(misrow(b), miscol(b)+1:end) = angle_source_correct(misrow(b), miscol(b)+1:end) - dtheta(misrow(b), miscol(b));
        angle_source_correct(misrow(b), miscol(b)+2:end) = angle_source_correct(misrow(b), miscol(b)+2:end) - dtheta(misrow(b), miscol(b)+1);
        dtheta(misrow(b), miscol(b)) = 0;
        dtheta(misrow(b), miscol(b)+1) = 0;
    end
end

if fig==1
    allerrs = unique(misrow);
    for f = 1:length(allerrs)
        disp([ num2str(length(misrow)) ' potential aberrant bouts found'])
        absc = miscol(misrow == allerrs(f));
        h = figure;
        subplot(1,3,1)
        plot(angle_source(allerrs(f),:))
        hold on
        plot(angle_source_correct(allerrs(f),:))
        plot(diff(angle_source(allerrs(f),:)) )
        plot(dtheta(allerrs(f),:))        
        plot(absc, angle_source_correct(allerrs(f),absc), '*')
        plot(dtheta(allerrs(f),:))
        title(num2str(allerrs(f)))
        subplot(1,3,2)
        plot3(1:length(x(1,:)), x(allerrs(f),:), y(allerrs(f),:))
        hold on
        plot3(absc, x(allerrs(f),absc), y(allerrs(f),absc), '*')
        hold off
        subplot(1,3,3)
        plot(x(allerrs(f),:), y(allerrs(f),:))
        hold on
        plot(x(allerrs(f),absc), y(allerrs(f),absc), '*')
        waitfor(h)
    end
end

end