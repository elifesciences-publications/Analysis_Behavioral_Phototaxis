function [X,Y,dX,dY ] = histprofile( data1, data2,numbin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dbin=(nanmax(data1)-nanmin(data1))/numbin;
X=(nanmin(data1)+dbin*0.5):dbin:(nanmax(data1)-dbin*0.5);
dX=dbin;
Y=[];
for n=1:numbin;
    w=find( data1>=(X(n)-dbin*0.5) & data1<(X(n)+dbin*0.5) );
    Y(n)=nanmean(data2(w));
    dY(n)=nanstd(data2(w))/sqrt(length(w));
end

end

