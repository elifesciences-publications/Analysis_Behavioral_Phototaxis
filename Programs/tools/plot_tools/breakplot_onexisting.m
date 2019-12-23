function h=BreakPlot_onexisting(x,y,y_break_start,y_break_end)
% BreakPlot(x,y,y_break_start,y_break_end,break_type)
% Produces a plot who's y-axis skips to avoid unnecessary blank space
%
% INPUT
% x
% y
% y_break_start
% y_break_end
% break_type
%    if break_type='RPatch' the plot will look torn
%       in the broken spaces
%    if break_type='Patch' the plot will have a more
%       regular, zig-zag tear
%    if break_plot='Line' the plot will merely have
%       some hash marks on the y-axis to denote the
%       break
%
% EXAMPLE #1:
%
% NEW = 1;
% sampleTimes = [1:200];
% RPM = [600+rand(1,100)*500, 4500+rand(1,100)*2500];
% figure;
% % BreakPlot
% str = {'','Line','Patch','RPatch'};
% for i=1:4
%    subplot(2,2,i);
%    if i==1
%        plot(sampleTimes,RPM,'r.')
%        title('Plot');
%    else
%        if ~NEW
%            breakplot_old(sampleTimes,RPM,2000,4000,str{i});
%        else
%            breakplot(sampleTimes,RPM,2000,4000,str{i},40);
%        end;
%        title(sprintf('Breakplot %s',str{i}));
%    end;
% end;
%
%
% EXAMPLE #2:
%
% figure;
% subplot(4,4,[1:2 5:6]);
% BreakPlot(rand(1,21),[1:10,40:50],10,40,'Line');
% subplot(4,4,[3:4 7:8]);
% BreakPlot(rand(1,21),[1:10,40:50],10,40,'Patch');
% subplot(4,4,[9:10 13:14]);
% BreakPlot(rand(1,21),[1:10,40:50],10,40,'RPatch');
% x=rand(1,21);y=[1:10,40:50];
% subplot(4,4,11:12);plot(x(y>=40),y(y>=40),'.');
% set(gca,'XTickLabel',[]);
% subplot(4,4,15:16);plot(x(y<=20),y(y<=20),'.');
%
%
% IT'S NOT FANCY, BUT IT WORKS.

% Michael Robbins
% michaelrobbinsusenet@yahoo.com

% TEST DATA
if nargin<4 y_break_end=39; end;
if nargin<3 y_break_start=11; end;
if nargin<2 y=[1:10,40:50]; end;
if nargin<1 x=rand(1,21); end;

% SOME DFINITIONS
y_break_width = y_break_end - y_break_start;
y_break_mid   = y_break_width./2 + y_break_start;
y_range       = range(y);

% LOSE THE DATA IN THE BREAK, WE DON'T NEED IT ANYMORE
i =  y>y_break_start & y <y_break_end;
x(i)=[];
y(i)=[];

% MAP THE DATA
i = y >= y_break_end;
y2 = y - i.*y_break_width;

% PLOT THE MAPPED DATA
% h    = plot(x,y,'-');
% h.LineWidth = 1.5;

h    = plot(x,y2,'-');
h.LineWidth = 1.5;

% CREATE THE "BREAK" EFFECT
xlim = get(gca,'xlim');
xtick      = get(gca,'XTick');
ytick      = get(gca,'YTick');
yticklabel = get(gca,'yticklabel');

y_gap_width = y_range ./ 100;
y_half_gap = y_gap_width./2;
y_gap_mid  = y_break_start + y_half_gap;


set(gca,'xlim',xlim);

% % MAP TICKS BACK
% i_wrong_ticks = ytick > y_break_start;
% ytick = ytick + i_wrong_ticks.*y_break_width;
% integer_ticks = all(floor(ytick) == ytick);
% label_width = size(yticklabel,2);
% if integer_ticks
%    format_string = sprintf('%%%dd\n',label_width);
% else
%    left_side = ceil(log10(max(ytick)));
%    right_side = label_width-left_side-1;
%    format_string = sprintf('%%%d.%df\n',label_width,right_side);
% end;
% set(gca, 'yticklabel', num2str(ytick'));