
% ----- Selection : on framerate (experiment removed if too slow) ---
if sum(D.experiment.framerate(:,3) < 20) > 1
    continue
end

% ----- Selection : on sequences length (short sequences removed)
[angleLab, angleSource, angleLabFilt , framerateInfo, coordinates, ~] = remove_sequence(D, timemin);
if isempty(angleSource)
    disp('to delete : ')
    disp([date fish])
    continue
end
framerate = framerateInfo(:,3);

% --- Correction : for potential aberrant bouts (head/tail misplacement)
fig = 0;
angleSourceCorr = find_aberrant_bouts(angleSource, coordinates, fig);

% --- Correct angle relative to source
if ~isfield(D.experiment, 'Illum')
    angleSourceCorr2 = angleSourceCorr + (angleLab(:,1) - angleSource(:,1));
    angleLabFilt = angleLabFilt + (angleLab(:,1) - angleSource(:,1));
else
    angleSourceCorr2 = angleSourceCorr;
    angleLabFilt = smoothCoord(angleSourceCorr2);
end

% --- Spot real bouts
fig = 1;
[boutIndices, ~, ~] = BoutSpot(coordinates, angleSourceCorr2, framerateInfo, fig);

del = input('sequences to delete manually : ');

if length(del) == size(angleSourceCorr,1)
    disp('whole experiment discarded : ')
    disp([fish date])
    continue
end
angleLab(del,:)=[];
angleSourceCorr2(del,:)=[];
angleLabFilt(del,:)=[];
framerateInfo(del,:)=[];
coordinates(:,:,del)=[];