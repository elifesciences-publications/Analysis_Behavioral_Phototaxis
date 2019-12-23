function [boutindexes] = BoutSpot(coordinates, angle_source, framerate)

coord = coordinates;
theta = angle_source;
fHz = framerate(:,3);
meanfhz = mean(fHz);
coord(coord==0) = NaN;
theta(theta ==0) = NaN; 

x = squeeze(coord(:,1,:))';
y = squeeze(coord(:,2,:))';

[thetaCorrect] = find_aberrant_bouts(theta, x, y, 1);

%smooth all variables
mtheta = NaN(size(thetaCorrect));
mx = NaN(size(x));
my = NaN(size(y));
for i = 1 : size(thetaCorrect,1)
    th = thetaCorrect(i,:);
    xf = x(i,:);
    yf = y(i,:);
    f = find(isnan(th),1);
    if isempty(f)
        mtheta(i,:) = smooth(th,7, 'sgolay')';
        mx(i,:) = smooth(xf, 11, 'sgolay')';
        my(i,:) = smooth(yf, 11, 'sgolay')';
    else
        mtheta(i,1:f-1) = smooth(th(1:f-1),7, 'sgolay');
        mx(i,1:f-1) = smooth(xf(1:f-1), 11, 'sgolay');
        my(i,1:f-1) = smooth(yf(1:f-1), 11, 'sgolay');
    end
    %plot(mtheta(i,1:f-1))
end

dx = diff(x, 1, 2);
dxcarr = dx.^2;

dy = diff(y, 1, 2);
dycarr = dy.^2;

dtheta = diff(thetaCorrect, 1,2);
dthetacarr = dtheta.^2;

vardth = nanvar(dtheta(:));
vardxy = nanvar(dx(:)+dy(:));

sigdisplacementmatrix = ((dthetacarr'/vardth).*((dxcarr'+dycarr')/vardxy))';
datasetSize = size(sigdisplacementmatrix,1);

%%
minIPI = 0.7; %sec
pkmaxnb = 0;
minh = 1;%prctile(sigdisplacementmatrix(:),80);
for i = 1 : datasetSize
    pks = findpeaks(sigdisplacementmatrix(i,:), 'MinPeakHeight', minh,...
        'MinPeakDistance', minIPI*fHz(i));
    if pkmaxnb < length(pks)
    pkmaxnb = length(pks);
    end
end
%%
maxbout = NaN(size(sigdisplacementmatrix,1),pkmaxnb);
boutduration = 0.4; %s
kin = NaN(size(sigdisplacementmatrix,1),pkmaxnb);

for i = 1 : datasetSize
    sdm = sigdisplacementmatrix(i,:);
    [pks, locs] = findpeaks(sdm, 'MinPeakHeight', minh,...
        'MinPeakDistance', minIPI*fHz(i));
    maxbout(i, 1:length(locs)) = locs;
    
    sdm(isnan(sdm))=0;
    se = strel('line', 7, 0);
    di = imdilate(sdm, se);
    er = imerode(di,se);
    
    btpts = round(boutduration*fHz(i));
    
    for j = 1 : length(locs)
%         if btpts < locs(j) && length(er) >= (locs(j) + btpts)
%             bout = er(locs(j) - btpts : locs(j) + btpts);
%             k = find(diff(bout)>0.02, 1);
%             if isempty(k)
%                 [~, k]=max(diff(bout));
%             end
%             kin(i,j) = k + locs(j)-1 - btpts;
%         elseif btpts < locs(j) && length(er) < (locs(j) + btpts)
%             bout = er(locs(j) - btpts : end);
%             [~, k] = find(diff(bout)>0.001, 1);
%             if isempty(k)
%                 [~, k]=max(diff(bout));
%             end
%             kin(i,j) = k + locs(j)-1 - btpts;
%         else
%             kin(i,j) = locs(j)-1;
%         end
%         test = 15; % pts
%         if kin(i,j) > test && length(sdm) - kin(i,j) > test
%             meanxb = nanmean(x(i,kin(i,j)-test:kin(i,j)));
%             meanyb = nanmean(y(i,kin(i,j)-test:kin(i,j)));
%             meanxa = nanmean(x(i,kin(i,j):kin(i,j)+test));
%             meanya = nanmean(y(i,kin(i,j):kin(i,j)+test));
%             meanthb = nanmean(mtheta(i,kin(i,j)-test:kin(i,j)));
%             meantha = nanmean(mtheta(i,kin(i,j):kin(i,j)+test));
%             if abs(meanxb-meanxa) < 2 && abs(meanyb-meanya) < 2 && abs(meanthb - meantha) < 3
%                 kin(i,j) = NaN;
%             end
%         end
    end %comment?
    
%     er(er==0) = NaN;
%     K = kin(i,1:length(locs));
%     replaceNaNs = length(K(isnan(K)));
%     K(isnan(K)) = [];
%     K = [K NaN(1,replaceNaNs)];
%     kin(i,1:length(locs)) = K;

    % plot
        figure
        plot(x(i,:))
        hold on
    plot(y(i,:))
    plot((theta(i,:)))
    plot((mtheta(i,:)))
    plot(sigdisplacementmatrix(i,:), 'k')
    plot(pks, locs, '*')
    plot(locs, x(i,locs),'o',locs, y(i,locs), 'o')
%     plot(kin(i,1:length(locs)-replaceNaNs), x(i,kin(i,1:length(locs)- replaceNaNs)), 'o')
%     plot(kin(i,1:length(locs)-replaceNaNs), y(i,kin(i,1:length(locs)- replaceNaNs)), 'o')
end

boutindexes = kin;
boutindexes(boutindexes==0) = NaN;
