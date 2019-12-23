
%1
ax = gca
x1 = ax.Children(2).XData;
y1 = ax.Children(2).YData;

[a1, b1] = cart2pol(x1, y1);
figure
a1 = unwrap(a1);
plot(diff(abs(a1))./cumsum(b1))

%%
c = imread('1.png');
im = 1-im2bw(c);
figure
imagesc(im)
colormap gray
axis image

figure
[n,r] = boxcount(im,'slope');

% The boxcount shows that the local exponent is approximately constant for
% less than one decade, in the range 8 < R < 128 (the 'exact' value of Df
% depends on the threshold, 80 gray levels here):

df = -diff(log(n))./diff(log(r));
disp(['Fractal dimension, Df = ' num2str(mean(df(4:8))) ' +/- ' num2str(std(df(4:8)))]);

hausDim(im)

