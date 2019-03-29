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
b=30;

[binvals_a, elts_per_bin_a, v2binMatrix_a] = BinsWithEqualNbofElements(Vart1_a, Vart2, b, b+10);
[binvals_b, elts_per_bin_b, v2binMatrix_b] = BinsWithEqualNbofElements(Vart1_b, Vart2, b, b+10);

%***
f = figure;
subplot(1,2,1)
errorbar(binvals_a,  mean(v2binMatrix_a, 2), std(v2binMatrix_a,1,2)/sqrt(elts_per_bin_a), 'k')
xlabel(xla)
ylabel(yl)
ax=gca;
ax.FontSize = 14;

subplot(1,2,2)
errorbar(binvals_b,  mean(v2binMatrix_b, 2), std(v2binMatrix_b,1,2)/sqrt(elts_per_bin_b), 'k')
xlabel(xlb)
ylabel(yl)
ax=gca;
ax.FontSize = 14;


%% Advancement
Vart1_a = dLum(:, 1:end-1);
Vart1_b = dLum(:, 1:end-1)./( (Lum(:, 1:end-2) + Lum(:, 2:end-1))./2 ) ;
xla = 'I';
xlb = 'dI/I';
Vart2 = Adv(:, 2:end)/pxmm;
b = 15;
yl = 'advancement';

[binvals_a, elts_per_bin_a, v2binMatrix_a] = BinsWithEqualNbofElements(Vart1_a, Vart2, b, b+10);
[binvals_b, elts_per_bin_b, v2binMatrix_b] = BinsWithEqualNbofElements(Vart1_b, Vart2, b, b+10);

%***
f = figure;
subplot(1,2,1)
errorbar(binvals_a,  median(v2binMatrix_a, 2), std(v2binMatrix_a,1,2)/sqrt(elts_per_bin_a), 'k')
xlabel(xla)
ylabel(yl)

subplot(1,2,2)
errorbar(binvals_b,  median(v2binMatrix_b, 2), std(v2binMatrix_b,1,2)/sqrt(elts_per_bin_b), 'k')
xlabel(xlb)
ylabel(yl)

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

subplot(1,2,2)
errorbar(binvals_b,  mean(v2binMatrix_b, 2), std(v2binMatrix_b,1,2)/sqrt(elts_per_bin_b), 'k')
xlabel(xlb)
ylabel(yl)

%% Velocity vs I & dI/I
Vart1_a = Lum(:, 1:end-2);
Vart1_b = dLum(:, 1:end-1);%./( (Lum(:, 1:end-2) + Lum(:, 2:end-1))./2 ) ;
xla = 'I';
xlb = 'dI/I';
Vart2 = D(:, 2:end)/pxmm./IBI(:, 2:end);
yl = 'Velocity (mm/s)';
b=30;

[binvals_a, elts_per_bin_a, v2binMatrix_a] = BinsWithEqualNbofElements(Vart1_a, Vart2, b, b+10);
[binvals_b, elts_per_bin_b, v2binMatrix_b] = BinsWithEqualNbofElements(Vart1_b, Vart2, b, b+10);

%***
f = figure;
subplot(1,2,1)
errorbar(binvals_a,  mean(v2binMatrix_a, 2), std(v2binMatrix_a,1,2)/sqrt(elts_per_bin_a), 'k')
xlabel(xla)
ylabel(yl)
ax=gca;
ax.FontSize = 14;


subplot(1,2,2)
errorbar(binvals_b,  mean(v2binMatrix_b, 2), std(v2binMatrix_b,1,2)/sqrt(elts_per_bin_b), 'k')
xlabel(xlb)
ylabel(yl)
ax=gca;
ax.FontSize = 14;