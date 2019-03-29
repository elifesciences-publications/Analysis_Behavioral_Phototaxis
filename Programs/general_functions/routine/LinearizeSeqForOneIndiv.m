function[SequencesLin] = LinearizeSeqForOneIndiv (M, rows)

separatedSequences = M(rows, :);
separatedSequencesCols = separatedSequences';
SequencesLin = separatedSequencesCols(:);
SequencesLin(isnan(SequencesLin)) = [];

end