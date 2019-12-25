%t=[1:10000];
dx=[];
psw_turn=0.1; % probabibility of switching left vs right states
p_turn=0.5; % probability of triggering a turn swim (otherwise go straight)
wturn=1;
wstraight=0.1;

%lrstate=(rand<pswitch)*2-1;
lrst=1;
dx=[];

 
for t=1:20000
    %pause(1)
    %lrstate=[lrstate,lrstate(end)*(-(rand<pswitch)*2+1)];
    lrst=lrst*((rand>psw_turn)*2-1);
    turn=rand<p_turn;
%     disp(['turn is ' num2str(turn)])
%     disp(['lrst is ' num2str(lrst)])
    if turn==1
        dxt=lrst*abs(wturn*randn);
    else
        dxt=lrst*abs(wstraight*randn);
    end
    dx=[dx,dxt];
end

 
plot(autocorr(dx),'*')