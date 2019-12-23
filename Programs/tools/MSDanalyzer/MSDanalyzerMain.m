% MSD analyzer

Xicell = cell(size(Xi,1),1);
for i = 1 : size(Xi,1)
    xsel = Xi(i,:);
    firstnan = find(isnan(xsel),1);
    if firstnan==1
        continue
    end
    mat = [find(xsel(1:firstnan-1))' xsel(1:firstnan-1)'];
    Xicell{i} = mat;
end
Xicell(cellfun(@isempty,Xicell)) = []

ma = msdanalyzer(1, 'rad', 'bouts');
ma = ma.addAll(Xicell);
ma = ma.computeMSD;
ma = ma.computeDrift('velocity');
ma.plotMSD;
ma.plotMeanMSD(gca, true)

