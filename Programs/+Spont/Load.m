
% .........................................................................
Es = Spont.load_spont_data;
% .........................................................................

% minimal bout number/fish threshold :
threshold = 30;
[Xi, FishNi, Ri, TimeBouti, xCoordi, yCoordi, FishID, sequencesperfish, boutsperfish]...
    = Spont.delete_f_not_enough_bouts(Es, threshold);

different_fish = unique(FishID);

% wrap & create dX
Xiwrapped = wrapToPi ( Xi ) ;
dXi = diff(Xi,1,2);
dXiw = diff(Xiwrapped, 1, 2);

% distance and projection on orientation vector
dxCoordi = diff(xCoordi, 1, 2);
dyCoordi = diff(yCoordi, 1, 2);
Disti = sqrt(dxCoordi.^2 + dyCoordi.^2);

% time
IBIi = diff(TimeBouti,1,2);

% trejctories
trajOrientationi = wrapToPi(atan(dyCoordi./dxCoordi));
trajOrientationi(dxCoordi < 0) = trajOrientationi(dxCoordi < 0) - pi;
Alpha = wrapToPi(Xi(:, 1:end-1)) + trajOrientationi;
Advi = Disti.*cos(Alpha);

% some stats per fish

fishdXstd = NaN(length(different_fish), 1);
fishdXmean = NaN(length(different_fish), 1);
fishBias = NaN(length(FishID),1);
fishDev = NaN(length(FishID),1);
scount = 1;
for i = different_fish'
    fishseqs = find(FishID == i);
    dxi = dXi(fishseqs,:);
    fishdXstd(i) = nanstd(dxi(:));
    fishdXmean(i) = nanmean(dxi(:));
    
    %--- vectors same size as dataset ---
    fishBias(scount : scount+length(fishseqs)-1) = repmat(fishdXmean(i), length(fishseqs), 1);
    fishDev(scount : scount+length(fishseqs)-1) = repmat(fishdXstd(i), length(fishseqs), 1);
    scount = scount + length(fishseqs);
end
dXissbiais = dXi - repmat(fishBias, [1 size(dXi,2)]);
dXin = dXissbiais ./ repmat(fishDev, [1 size(dXi,2)]);
dXissdev = dXi ./ repmat(fishDev, [1 size(dXi,2)]);