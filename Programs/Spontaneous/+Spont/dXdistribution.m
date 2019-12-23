<<<<<<< HEAD
function[Wturn, Wfor, Pturn] = dXdistribution(dXi, different_fish)
=======
function[Wturn, Wfor, Pturn] = dXdistribution(dXi, FishID, different_fish)
>>>>>>> master

% Compare built-in and custom fits

[colour] = colour_palette(0, 1);

<<<<<<< HEAD
%% Prepare data (into bins)

V = dXi(:);
V(isnan(V))=[];
N = length(V);
binwidth = 3*iqr(dXi(:))/(N^(1/3));
bins = min(dXi(:)) : binwidth : max(dXi(:));
[yhist, edges]  = histcounts(V, bins);
xhist = edges(1:end-1) + binwidth/2;
yhist = yhist/(sum(yhist)*binwidth);

%% Fits and Plots
% --- built-in gauss 2 on whole dataset ---
f = fit(xhist.',yhist.','gauss2', 'Robust', 'Bisquare');

%***
fig1 = figure;
fig1.Name = 'built-in fit';
semilogy(xhist, yhist)
hold on
semilogy(xhist, f.a1*exp(-((xhist-f.b1)/f.c1).^2) + f.a2*exp(-((xhist-f.b2)/f.c2).^2))

xlim([-2*pi/3, 2*pi/3])
ylim([10e-5 10])
xticks([-pi/2,-pi/4, 0,pi/4, pi/2])
xticklabels({'-\pi/2','-\pi/4', '0','\pi/4', '\pi/2'})

=======
%%
% --- built-in gauss 2 on whole dataset ---
N = dXi(:);
N(isnan(N))=[];
N = length(N);
binwidth = 3*iqr(dXi(:))/(N^(1/3));
bins = min(dXi(:)):binwidth:max(dXi(:));
[yhist, edges]  = histcounts(dXi(:),bins);
xhist=edges(1:end-1)+binwidth/2;
yhist=yhist/(sum(yhist)*binwidth);

f = fit(xhist.',yhist.','gauss2');
forward = @(x) f.a1*exp(-((x-f.b1)/f.c1).^2);
side = @(x) f.a2*exp(-((x-f.b2)/f.c2).^2);

%***
figure;
semilogy(xhist, yhist, 'k', 'LineWidth', 1.5)
hold on
semilogy(xhist, f.a1*exp(-((xhist-f.b1)/f.c1).^2) + f.a2*exp(-((xhist-f.b2)/f.c2).^2),...
     'LineWidth', 1.5, 'Color', colour(5,:))
xlim([-2*pi/3, 2*pi/3])
ylim([10e-5 10])
>>>>>>> master
t = text(-pi/6, 0.01, {['\sigma_t = ' num2str(f.c2, 2) ' rad'];...
    ['\sigma_f = ' num2str(f.c1,2) ' rad'];...
    ['p_t = ' num2str(1-f.a1*f.c1/(f.a1*f.c1+f.a2*f.c2),2)];...
    ['p_f = ' num2str(f.a1*f.c1/(f.a1*f.c1+f.a2*f.c2))]}) ;
t.FontSize = 14;
t.FontName = 'Times New Roman';
title('from built-in fit')

<<<<<<< HEAD
legend(['data : ' num2str(different_fish(end)) ' fish'], 'N (\mu_t+\mu_f, \sigma_t+\sigma_f)')

ax = gca;
=======
%%
% --- custom gauss 2 ---
[f, xhist, yhist, muabs, w, Wturn, Wfor, Pturn] = Gauss2custom_simple(dXi);
wturn = Wturn(1);
wfor = Wfor(1);
pturn = Pturn(1);

%***
fig = figure;
plot(f, xhist',yhist');
ax = gca;

>>>>>>> master
ax.YScale = 'log';
ax.LineWidth = 2;
ax.FontSize = 16;
ax.FontName = 'Times New Roman';
ax.Children(1).Color = colour(4,:);
ax.Children(1).LineWidth = 2;
ax.Children(2).Color = [0 0 0];
ax.Children(2).Marker = 'none';
ax.Children(2).LineStyle = '-';
ax.Children(2).LineWidth = 1.5;
<<<<<<< HEAD
ax.Children(3).Color = colour(4,:);
ax.Children(3).Marker = 'none';
ax.Children(3).LineStyle = '-';
ax.Children(3).LineWidth = 1.5;
ax.TickLength = [0 0];

%%
% --- custom gauss 2 ---
[fc, xhist, yhist, muabs, w, Wturn, Wfor, Pturn] = Gauss2custom_simple(V);
wturn = Wturn(1);
wfor = Wfor(1);
pturn = Pturn(1);

%***
fig2 = figure;
fig2.Name = 'custom fit';
plot(fc, xhist',yhist');

=======
>>>>>>> master
ylim([10e-5 10])
xlim([-2*pi/3, 2*pi/3])
xticks([-pi/2,-pi/4, 0,pi/4, pi/2])
xticklabels({'-\pi/2','-\pi/4', '0','\pi/4', '\pi/2'})

xlabel('\delta\theta')
ylabel('pdf')
    
hold on
<<<<<<< HEAD
t = text(-pi/6, 0.01, {['\sigma_t = ' num2str(Wturn(1), 2)];...
    ['\sigma_f = ' num2str(Wfor(1),2) ' rad'];...
    ['p_t = ' num2str(Pturn(1),2)];...
    ['p_f = ' num2str(1-Pturn(1),2) ]}) ;
t.FontSize = 14;
t.FontName = 'Times New Roman';
title('from custom fit')

legend(['data : ' num2str(different_fish(end)) ' fish'], 'N (0, \sigma_t+\sigma_f)')

ax = gca;
ax.YScale = 'log';
ax.LineWidth = 2;
ax.FontSize = 16;
ax.FontName = 'Times New Roman';
ax.Children(1).Color = colour(4,:);
ax.Children(1).LineWidth = 2;
ax.Children(2).Color = [0 0 0];
ax.Children(2).Marker = 'none';
ax.Children(2).LineStyle = '-';
ax.Children(2).LineWidth = 1.5;
ax.Children(3).Color = colour(4,:);
ax.Children(3).Marker = 'none';
ax.Children(3).LineStyle = '-';
ax.Children(3).LineWidth = 1.5;
ax.TickLength = [0 0];

%%
muf = 0;
myfittype = fittype('sum_gauss_constr(x, a, mut, sigmat, muf, sigmaf)',...
    'coefficients', {'a','mut','sigmat', 'sigmaf'}, 'problem', {'muf'});
myfit = fit(xhist',yhist',myfittype, 'startpoint', [0.5 0 0.7 0.1], 'Upper', [1 Inf 0.9 1], 'problem', {muf});

%***
figure
plot(myfit, xhist, yhist)
=======
t = text(-pi/6, 0.01, {['\sigma_t = ' num2str(Wturn(1), 2) ' [ \pm' num2str(Wturn(2)-Wturn(1),2) ']_{ci99}' ' rad'];...
    ['\sigma_f = ' num2str(Wfor(1),2) ' [ \pm' num2str(Wfor(2)-Wfor(1),2) ']_{ci99}' ' rad'];...
    ['p_t = ' num2str(Pturn(1),2) ' [ \pm' num2str(Pturn(1)-Pturn(2),2) ']_{ci99}'];...
    ['p_f = ' num2str(1-Pturn(1),2) ' [ \pm' num2str(Pturn(1)-Pturn(2),2) ']_{ci99}']}) ;
t.FontSize = 14;
t.FontName = 'Times New Roman';
title('from custom-in fit')

legend(['data : ' num2str(different_fish(end)) ' fish'], 'N (\mu_t +\mu_f, \sigma_t+\sigma_f)')

%%
v = (dXi(:));
v = v(:);
v(isnan(v)) = [];
N = length(v);
binwidth = 3*iqr(v)/(N^(1/3));
bins = min(v) : binwidth : max(v);
[dx_histogram, x_histogram]  = histcounts(v,bins);
x_histogram = x_histogram(1:end-1) + binwidth/2;
dx_histogram = dx_histogram/(sum(dx_histogram)*binwidth);
muf = 0;
sigmaf = 0.095;
myfittype = fittype('sum_gauss_constr(x, a, mut, sigmat, muf, sigmaf)',...
    'coefficients', {'a','mut','sigmat'}, 'problem', {'muf', 'sigmaf'}); %
myfit = fit(x_histogram',dx_histogram',myfittype, 'startpoint', [0.5 0 0.7], 'Upper', [1 Inf 0.9], 'problem', {muf, sigmaf});%

%***
figure
plot(myfit, x_histogram, dx_histogram)
>>>>>>> master
myfit

