function [mtheta, mx, my] = smoothCoord(theta, x, y)

switch nargin
    case 1
        mtheta = NaN(size(theta));
        for i = 1 : size(theta,1)
            th = theta(i,:);
            f = find(isnan(th),1);
            if isempty(f)
                mtheta(i,:) = smooth(th, 5, 'sgolay')';
            else
                mtheta(i,1:f-1) = smooth(th(1:f-1),5, 'sgolay');
            end
        end
        mx = [];
        my = [];
    case 3
        mtheta = NaN(size(theta));
        mx = NaN(size(x));
        my = NaN(size(y));
        for i = 1 : size(theta,1)
            th = theta(i,:);
            xf = x(i,:);
            yf = y(i,:);
            f = find(isnan(th),1);
            if isempty(f)
                mtheta(i,:) = smooth(th, 5, 'sgolay')';
                mx(i,:) = smooth(xf, 12, 'sgolay')';
                my(i,:) = smooth(yf, 12, 'sgolay')';
            else
                mtheta(i,1:f-1) = smooth(th(1:f-1),5, 'sgolay');
                mx(i,1:f-1) = smooth(xf(1:f-1), 12, 'sgolay');
                my(i,1:f-1) = smooth(yf(1:f-1), 12, 'sgolay');
            end
        end
    case 2
        disp('in smoothCorrd : x or y missing')
end
