
% All Rs

%p1 = sin6

R = {'lat', 0.2120, 0.1908;
    'p1', 0.1034, 0.1198 ;
    'p2', 0.2017, 0.0996 ;
    'p3', 0.1371, 0.0874 ;};

Rmat = cell2mat(R(:,2:end));

figure
bar(Rmat)
xticklabels({'lat' 'p1' 'p2' 'p3'})
legend('data', 'simulation')
ylabel('R.cos(<\theta>)')

ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;