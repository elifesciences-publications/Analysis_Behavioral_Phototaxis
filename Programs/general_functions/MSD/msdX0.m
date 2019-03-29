function [msd0] = msdX0(X, bias)

% INPUT X
%   rows : sequences
%   columns : time
if size(bias,1) == 1 && bias ~= 0 
    b = 0 : bias : bias*(size(X,2)-1);
    b =  repmat(b,size(X,1),1);
elseif size(bias,1) == size(X,1)
    b = NaN(size(X));
    for i = 1 : size(X,2)
        b(:,i) = bias*(i-1);
    end
else
    b = 0;
end

X_X0 = X - X(:,1) - b;
sd = (X_X0).^2;
msd0 = mean(sd, 1, 'omitnan');

end