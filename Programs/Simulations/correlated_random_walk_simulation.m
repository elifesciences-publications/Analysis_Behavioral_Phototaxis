
%% random walk simulations
rep = 1000;
VAR = NaN(rep,1);
a = NaN(rep,1);
intac = NaN(rep,1);
stdev = NaN(rep,1);

for i = 1 : rep
    
    % distribution size parameters
    wturn = 0.59;
    wfor = 0.093;
    pturn = 0.4143;
    psw_turn = rand;
    
    %%
    rng('shuffle')
    Nexp = 1000; % number of experiments
    Ntimes = 100; % number of time steps per experiment
    theta_complete = NaN(Nexp, Ntimes); % complete dataset
    
    theta_ini = rand(1,Nexp)*2*pi; % initial angles uniformly distributed
    
    lrst = (rand > 0.5)*2-1; % intial turn state random 50/50
    
    for N = 1:Nexp
        theta = NaN(1,Ntimes);
        theta(1) = theta_ini(N);
        for t = 1:Ntimes-1
            th = theta(t);
            lrst = lrst*((rand > psw_turn)*2-1);
            turn = rand < pturn;
            if turn
                dth = lrst * abs(wturn*randn);
            else
                dth = lrst * abs(wfor*randn);
            end
            theta(t+1) = th + dth;
        end
        theta_complete(N,:) = theta;
    end
    
    dth_complete = diff(theta_complete,1,2);
    
    %% test D measure
    VAR(i) = psw_turn;
    msdth = msdX0shuffled(theta_complete,0);
    lf = polyfit((10:30), msdth(10:30),1);
    a(i) = lf(1);
  
    ac = xcorrMatrixRows (dth_complete);
    intac(i) = sum(ac);
    stdev(i) = std(dth_complete(:));
    disp(i)
end

%***
figure
plot([0.1 1], [0.1 1], 'k--', 'DisplayName', 'y=x')
hold on
ps = 0:0.05:1;
cm = colormap(jet);
for j = 1:length(ps)-1
    pswvalues = find(VAR>ps(j) & VAR < ps(j+1));
    color = cm(round(length(cm)*j/length(ps)),:);
    if ps(j) == 0 || round(ps(j)*1000)/1000 == 0.15 || ps(j) == 0.5
        disp('y')
        plot(a(pswvalues), intac(pswvalues).*(stdev(pswvalues)).^2, 'o',...
            'MarkerFaceColor', color, 'MarkerEdgeColor', color,...
            'DisplayName', [num2str(ps(j)) '-' num2str(ps((j+1)))])
    else
        plot(a(pswvalues), intac(pswvalues).*(stdev(pswvalues)).^2, 'o',...
            'MarkerFaceColor', color, 'MarkerEdgeColor', color,...
            'HandleVisibility','off')
    end
end
ax=gca;
ax.XScale = 'log';
ax.YScale = 'log';
grid on
legend
xlabel('a from ax+b fit on MSD')
ylabel('var(dx).\int acf')
title('different Pswitch values')

y= (intac.*(stdev).^2) ./ a;
x = VAR;
ft = fittype('exp(a*x)+exp(b*sqrt(x))+c', 'coefficients', {'a','b','c'});
myfit = fit(x,y,ft, 'StartPoint', [-60 1.2 -0.8]);

%***
figure
plot(myfit, x, y)
xlabel('P switch')
ylabel('<\delta\theta^2> \int acf / a')
legend('data from simu', 'fitted curve')

root_path = '/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Analysis/';
folder = 'SimuData';
SimuPswitch.fit = myfit;
SimuPswitch.x = x;
SimuPswitch.y = y;
mkdir([root_path folder])
save([root_path folder filesep 'SimuPswitch.mat'], 'SimuPswitch')



