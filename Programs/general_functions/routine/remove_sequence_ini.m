function [angle, coordinates, framerate, seq_remove]=remove_sequence_ini(D,timemin)
%% Comments
% Select the sequences to analyse. Remove the empty sequences and the
% sequence shorter than time
%Inputs -----
% D: structure from the experiment
% time: time minimal of the sequence for analysing, in seconds
%Outputs -----
% seq_remove: list of the sequences which have been removed
% sequence: list of the remaining sequences
% framerate, coordinates, angle: data with only the remaining sequences

% /!\ Specific for the initilization analysis

%% Code
if isfield(D, 'initialization')
    angle = D.initialization.angleIni;
    framerate = D.initialization.framerateIni;
    coordinates = D.initialization.coordinatesIni;
else
    angle = D.angleIni;
    framerate = D.framerateIni;
    coordinates = D.coordinatesIni;
end

seq_remove = [];
sequence = linspace(1,size(angle,1),size(angle,1));
for seq = 1:size(angle,1)
    % ----- remove empty sequence -----
    a = nansum(angle(seq,:));
    if a == 0
        seq_remove = [seq_remove seq];
    else
        endseq = framerate(seq,4);
        if endseq == 0
            continue
        end
        if endseq > size(angle,2)
            endseq = size(angle,2);
        end
        angle(seq,endseq:end) = nan;
        coordinates(seq,endseq:end) = nan;
    end
    % ----- remove sequence shorter than time second -----
    t = framerate(seq,2) - framerate(seq,1);
    if t < timemin
        seq_remove = [seq_remove seq];
    end  
end

seq_remove = unique(seq_remove);
sequence = setdiff(sequence,seq_remove);

angle(seq_remove,:) = [];
framerate(seq_remove,:) = [];
coordinates(:,:,seq_remove) = [];