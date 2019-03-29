function [MSD,stdMSD,n] = msdX0shuffled2( X, bias )

[n, m] = size( X );

if size(bias,1) == 1 && bias ~= 0 
    b = 0 : bias : bias*(m-1);
    b =  repmat(b,size(X,1),1);
elseif size(bias,1) == n
    b = NaN(n,m);
    for i = 1 : m
        b(:,i) = bias*(i-1);
    end
else
    warning('bias set to 0')
    b = 0;
end

squared_displacement = NaN(n*m,m);
X_Xt0 = X - X(:,1) - b;
squared_displacement1 = (X_Xt0).^2;
squared_displacement(1:n,1:m) = squared_displacement1;

for j = 2 : m
    XiShift = X(:,j:end);
    if size(bias,1) == 1 && bias ~= 0 
        b = 0 : bias : bias*(size(XiShift,2)-1);
        b =  repmat(b,size(XiShift,1),1);
    elseif size(bias,1) == size(XiShift,1)
        b = NaN(size(XiShift));
        for i = 1 : size(XiShift,2)
            b(:,i) = bias*(i-1);
        end
    else
        b = 0;
    end
    XiShift_Xt0 = XiShift - XiShift(:,1) - b;
    msd = ( XiShift_Xt0 ).^2;
    squared_displacement( (j-1) * n+1 : j*n , 1 : m-j+1 ) = msd;
end

deleteNaNrows = sum( isnan( squared_displacement ), 2 );
squared_displacement( deleteNaNrows == m, : ) = [];
MSD = mean(squared_displacement, 1, 'omitnan');
stdMSD = nanstd(squared_displacement,1,1);
n = size(squared_displacement,1)-sum(isnan(squared_displacement),1);

end