function[Xi, TimeBouti, xCoordi, yCoordi, FishNi] = all_experiments_pool(timemin, minBoutNumber)

% initial LOOP
% find the approximate number of bouts to be spotted per sequence
% (overestimation)

N=4;

L = 0 ;
l = 0 ;
Larg = NaN(N,1);
long = NaN(N,1);
for EXP = 1:N
    if EXP == 1
        [Dates, Fish, exptype, ~, intmax] = Temp.perexperiment.choose_experiment('exp60');
    elseif EXP == 2
        [Dates, Fish, exptype, ~, intmax] = Temp.perexperiment.choose_experiment('exp30');
    elseif EXP ==3
        [Dates, Fish, exptype, ~, intmax] = Temp.perexperiment.choose_experiment('sin60');
    elseif EXP == 4
        Lat.DatesFish;
        intmax = 0;
    end
   
    [Larg(EXP),long(EXP)] = Spont.spont_initial_loop(Dates, Fish, exptype, intmax, timemin);
    L = L+Larg(EXP);
    l = l+long(EXP);
    disp(Dates)
    disp(l)
end

%% POOLING LOOP
Xi = NaN(l,L);
TimeBouti = NaN(l,L);
xCoordi = NaN(l,L);
yCoordi = NaN(l,L);
FishNi = NaN(l,1);

row = 1;
for EXP = 1:N
    if EXP == 1
        [Dates, Fish, exptype, ~, intmax] = Temp.chooseExpType('exp60');
    elseif EXP == 2
        [Dates, Fish, exptype, ~, intmax] = Temp.chooseExpType('exp30');
    elseif EXP == 3
        [Dates, Fish, exptype, ~, intmax] = Temp.chooseExpType('sin60');
    elseif EXP == 4
        Lat.DatesFish
        intmax = 0;
    end
    
    [Xi_e, TimeBouti_e, xCoordi_e, yCoordi_e, FishNi_e] = ...
    Spont.spont_pooling_loop(Dates, Fish, long(EXP), Larg(EXP), exptype, intmax, timemin, minBoutNumber);
    
    sizeExp = size(Xi_e);
        
    Xi( row : row + sizeExp(1)-1 , 1 : sizeExp(2)) = Xi_e;
    TimeBouti( row : row + sizeExp(1)-1 , 1 : sizeExp(2)) = TimeBouti_e;
    xCoordi( row : row + sizeExp(1)-1 , 1 : sizeExp(2)) = xCoordi_e;
    yCoordi( row : row + sizeExp(1)-1 , 1 : sizeExp(2)) = yCoordi_e;
    FishNi( row : row + sizeExp(1)-1, 1 ) = FishNi_e;
    
    row = row +  sizeExp(1);
end

% clean up & tranform for useful variables
dxCoordi = diff(xCoordi, 1, 2);
dyCoord = diff(yCoordi, 1, 2);
Disti = sqrt(dxCoordi.^2+ dyCoord.^2);

% calculate advancement
trajOrientation = wrapToPi(atan(dyCoord./dxCoordi));
trajOrientation(dxCoordi < 0) = trajOrientation(dxCoordi < 0) - pi;
dAlpha = wrapToPi(Xi(:, 1:end-1)) + trajOrientation;

Ri = Disti.*cos(dAlpha);

% completely delete sequences where cumulative advancement <0
cumulativeAdvancement = nansum(Ri, 2);
cumathresh = 100;
seqs_to_delete = find(cumulativeAdvancement < cumathresh);

Xi(seqs_to_delete, :) = []; 
TimeBouti(seqs_to_delete, :) = []; 
xCoordi(seqs_to_delete, :) = []; 
yCoordi(seqs_to_delete, :) = []; 
FishNi(seqs_to_delete, :) = []; 

% replace aberrant points by NaNs
dXi = diff( Xi, 1, 2 );

dxCoordi = diff(xCoordi, 1, 2);
dyCoord = diff(yCoordi, 1, 2);
Disti = sqrt(dxCoordi.^2+ dyCoord.^2);

% calculate advancement
trajOrientation = wrapToPi(atan(dyCoord./dxCoordi));
trajOrientation(dxCoordi < 0) = trajOrientation(dxCoordi < 0) - pi;
dAlpha = wrapToPi(Xi(:, 1:end-1)) + trajOrientation;

Ri = Disti.*cos(dAlpha);
Ti = Disti.*sin(dAlpha);

% delete aberrant points
[S, B] = find(Ri<0 & abs(dXi)<1);
Xi(S,B) = NaN;
TimeBouti(S,B) = NaN;
xCoordi(S,B) = NaN;
yCoordi(S,B) = NaN;

% --- orientation ---
Xi(Xi==0) = NaN;
[Xi] = deleteEmptyRows(Xi);
[TimeBouti] = deleteEmptyRows(TimeBouti);
[xCoordi] = deleteEmptyRows(xCoordi);
[yCoordi] = deleteEmptyRows(yCoordi);
[FishNi] = deleteEmptyRows(FishNi);

