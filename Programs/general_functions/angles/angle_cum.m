function [angle] = angle_cum(anglewrapped, framerate)
%% Comments
% Remove the 0-360? edges and smooth the angle
%Inputs -----
% angleini: raw angle from the initialization
% framerate: framerate of the initialization
% fig: if fig=0, don't plot, if fig=1, plot
%Outputs -----
% angle: raw angle with no 0-360 edge
% angle_mov: movmean of angle
% ang_rand: randomize angle with no 0-360 edge (translation of the raw
% data, not the same angle(0))
% mov_ang_rand: movmean of ang_rand


%% Code
angle = nan(size(anglewrapped));

for seq = 1:size(anglewrapped,1)
    endseq = framerate(seq,4);
    if endseq > size(anglewrapped,2)
        endseq = size(anglewrapped,2);
    end
    angle(seq,1) = anglewrapped(seq,1);
    for j = 2:endseq-1
        d = anglewrapped(seq,j)-anglewrapped(seq,j-1);
        d = angle_per_frame(d);
        if isnan(d)
            continue
        elseif abs(d) < 150
            angle(seq,j) = angle(seq,j-1) + d;
        else
            angle(seq,j) = angle(seq,j-1);
        end
    end
end