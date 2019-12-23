function[col] = sizeWithoutNans(M)

dim = 2; % 2: first NaN in row

[~, col] = max( isnan(M), [], dim );

end