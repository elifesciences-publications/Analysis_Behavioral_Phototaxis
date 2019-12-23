
% All Rs

% bouts 2 to med
R = {'lat', 0.2120, 0.1908;
    'unif 1', 0.1034, 0.1298 ;
    'unif 2', 0.2017, 0.1364 ;
    'unif 3', 0.1371, 0.1164 ;};

% % bouts 8 to med
% R = {'lat', 0.2120, 0.1908;
%     'unif 1', 0.1034, 0.1915 ;
%     'unif 2', 0.2017, 0.1236 ;
%     'unif 3', 0.128, 0.0869 ;};

Rmat = cell2mat(R(:,2:end));

figure
bar(Rmat)
xticklabels({'lat' 'p1' 'p2' 'p3'})
legend('data', 'simulation')
ylabel('R.cos(<\theta>)')

ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;