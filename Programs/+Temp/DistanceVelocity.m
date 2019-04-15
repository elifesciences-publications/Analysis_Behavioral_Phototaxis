%% Distances and velocities between bouts in temporal phototaxis

% Distance : distance between two bouts

% Advancement : distance between two bouts projected on the orientation
%               vector of the fish at previous bout

% T : transversal displacement, distance projected on vector orthogonal to
%     the fish's orrientation at previous bout

% Velocity : not really the instantaneous velocity, bout rather the distance between two
%           bouts divided by the time between two bouts

%% Distance vs I & dI/I
Vart1_a = Lum(:, 1:end-2);
Vart1_b = dLum(:, 1:end-1);%./( (Lum(:, 1:end-2) + Lum(:, 2:end-1))./2 ) ;
xla = 'I';
xlb = 'dI/I';
Vart2 = D(:, 2:end)/pxmm;
yl = 'IB distance (mm)';
b=7;

[binvals_a, elts_per_bin_a, v2binMatrix_a] = BinsWithEqualNbofElements(Vart1_a, Vart2, b, b+10);
[binvals_b, elts_per_bin_b, v2binMatrix_b] = BinsWithEqualNbofElements(Vart1_b, Vart2, b, b+10);

%***
f = figure;
subplot(1,2,1)
errorbar(binvals_a,  mean(v2binMatrix_a, 2), std(v2binMatrix_a,1,2)/sqrt(elts_per_bin_a), 'k')
xlim([0 0.25])
ylim([2.3 3.3])
xlabel(xla)
ylabel(yl)
ax=gca;
ax.FontSize = 18;

subplot(1,2,2)
errorbar(binvals_b,  mean(v2binMatrix_b, 2), std(v2binMatrix_b,1,2)/sqrt(elts_per_bin_b), 'k')
xlim([-0.1 0.1])
ylim([2.3 3.3])
xlabel(xlb)
ylabel(yl)
ax=gca;
ax.FontSize = 18;


%% Advancement
Vart1_a = Lum(:, 1:end-2);
Vart1_b = dLum(:, 1:end-1);
Vart1_c = dLum(:, 1:end-1)./( (Lum(:, 1:end-2) + Lum(:, 2:end-1))./2 ) ;
xla = 'I';
xlb = 'dI';
xlc = 'dI/I';
Vart2 = Adv(:, 2:end)/pxmm;
b = 12;
yl = 'advancement';

[binvals_a, elts_per_bin_a, v2binMatrix_a] = BinsWithEqualNbofElements(Vart1_a, Vart2, b, b+3);
[binvals_b, elts_per_bin_b, v2binMatrix_b] = BinsWithEqualNbofElements(Vart1_b, Vart2, b, b+3);
[binvals_c, elts_per_bin_c, v2binMatrix_c] = BinsWithEqualNbofElements(Vart1_c, Vart2, b, b+3);


%***
f = figure;
subplot(1,3,1)
errorbar(binvals_a,  median(v2binMatrix_a, 2), std(v2binMatrix_a,1,2)/sqrt(elts_per_bin_a), 'k', 'LineWidth', 1.5)
ylim([1.6 2.8])
xlabel(xla)
ylabel(yl)
ax=gca;
ax.FontSize = 18;

subplot(1,3,2)
errorbar(binvals_b,  median(v2binMatrix_b, 2), std(v2binMatrix_b,1,2)/sqrt(elts_per_bin_b), 'k', 'LineWidth', 1.5)
ylim([1.6 2.8])
xlabel(xlb)
ylabel(yl)
ax=gca;
ax.FontSize = 18;

subplot(1,3,3)
errorbar(binvals_c,  median(v2binMatrix_c, 2), std(v2binMatrix_c,1,2)/sqrt(elts_per_bin_c), 'k', 'LineWidth', 1.5)
ylim([1.6 2.8])
xlabel(xlc)
ylabel(yl)
ax=gca;
ax.FontSize = 18;

%% T vs I & dI/I
Vart1_a = Lum(:, 1:end-2);
Vart1_b = dLum(:, 1:end-1)./( (Lum(:, 1:end-2) + Lum(:, 2:end-1))./2 ) ;
xla = 'I';
xlb = 'dI/I';
Vart2 = T(:, 2:end)/pxmm;
b = 15;
yl = 'lateral displacement';

[binvals_a, elts_per_bin_a, v2binMatrix_a] = BinsWithEqualNbofElements(Vart1_a, Vart2, b, b+10);
[binvals_b, elts_per_bin_b, v2binMatrix_b] = BinsWithEqualNbofElements(Vart1_b, Vart2, b, b+10);

%***
f = figure;
subplot(1,2,1)
errorbar(binvals_a,  mean(v2binMatrix_a, 2), std(v2binMatrix_a,1,2)/sqrt(elts_per_bin_a), 'k')
xlabel(xla)
ylabel(yl)
ax=gca;
ax.FontSize = 18;

subplot(1,2,2)
errorbar(binvals_b,  mean(v2binMatrix_b, 2), std(v2binMatrix_b,1,2)/sqrt(elts_per_bin_b), 'k')
xlabel(xlb)
ylabel(yl)
ax=gca;
ax.FontSize = 18;

%% Velocity vs I & dI/I
Vart1_a = Lum(:, 1:end-2);
Vart1_b = dLum(:, 1:end-1);%./( (Lum(:, 1:end-2) + Lum(:, 2:end-1))./2 ) ;
Vart1_c = dLum(:, 1:end-1)./( (Lum(:, 1:end-2) + Lum(:, 2:end-1))./2 ) ;
xla = 'I';
xlb = 'dI';
xlc = 'dI/I';
Vart2 = D(:, 2:end)/pxmm./IBI(:, 2:end);
yl = 'Mean velocity (mm/s)';
b=12;

[binvals_a, elts_per_bin_a, v2binMatrix_a] = BinsWithEqualNbofElements(Vart1_a, Vart2, b, b+3);
[binvals_b, elts_per_bin_b, v2binMatrix_b] = BinsWithEqualNbofElements(Vart1_b, Vart2, b, b+3);
[binvals_c, elts_per_bin_c, v2binMatrix_c] = BinsWithEqualNbofElements(Vart1_c, Vart2, b, b+3);

%***
f = figure;
subplot(1,3,1)
errorbar(binvals_a,  mean(v2binMatrix_a, 2), std(v2binMatrix_a,1,2)/sqrt(elts_per_bin_a), 'k', 'LineWidth', 1.5)
xlabel(xla)
ylabel(yl)
ax=gca;
ax.FontSize = 18;

subplot(1,3,2)
errorbar(binvals_b,  mean(v2binMatrix_b, 2), std(v2binMatrix_b,1,2)/sqrt(elts_per_bin_b), 'k', 'LineWidth', 1.5)
xlabel(xlb)
ylabel(yl)
xlim([-0.1 0.1])
ax=gca;
ax.FontSize = 18;

subplot(1,3,3)
errorbar(binvals_c,  mean(v2binMatrix_c, 2), std(v2binMatrix_c,1,2)/sqrt(elts_per_bin_c), 'k', 'LineWidth', 1.5)
xlabel(xlc)
ylabel(yl)
ax=gca;
ax.FontSize = 18;