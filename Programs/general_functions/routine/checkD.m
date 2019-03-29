
% check D structure
X = [];
f = [];
count = 1
for p = 1: length(Dates)
    for q = 1: length(Fish(p,:))
        if iscellstr(Fish(p,q)) == 0
            continue
        end
        fish = char(Fish(p,q));
        date = char(Dates(p));
        
        
        load(['/Users/karp/Documents/PhD/Projects/Behaviorfish/PhototaxisFreeSwim/Data/'...
            exptype '/' date '/' fish '/data' fish '.mat']);
        
        if isfield(D.experiment, 'Illum')
            x = D.experiment.angleCum;
            x(x==0) = NaN;
            X = [X ; x(:)];
            f= [f, count];
            disp([date fish])
        end
        count = count+1;
    end
end 
