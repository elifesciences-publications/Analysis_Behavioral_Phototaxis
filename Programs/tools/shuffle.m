function B = shuffle(A)
randomIdxsA = randperm(numel(A));
B = reshape(A(randomIdxsA),size(A));
end