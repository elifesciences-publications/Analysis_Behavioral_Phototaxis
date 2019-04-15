<<<<<<< HEAD
function[] = auto_co_reinforcement(dXl,dXissbiais, DIlr, bins1, bins2, forceonzero, contrast)

% --- appearance ---
=======
function[] = auto_co_reinforcement(dXl,dXissbiais, DIlr, bins1, bins2, forceonzero)

>>>>>>> master
colour = colour_palette(0,3);
col_reinf = colour(3,:);
col_inhib = colour(4,:);

<<<<<<< HEAD
% --- 
=======
>>>>>>> master
Cnp1 = DIlr(:, 2:end-1)/max(abs(DIlr(:)));
dXn = dXl(:,1:end-1);
dXnp1 = dXl(:,2:end);
Prodn_np1 = (dXn.*dXnp1);
Absn_np1 = (abs(dXn).*abs(dXnp1));
ACn_np1 = Prodn_np1./Absn_np1;

limI = 0;
limX = 0;
reinf = find( (Cnp1 > limI & dXn >limX) | (Cnp1 < limI & dXn <-limX) );
inhib = find( (Cnp1 < limI & dXn >limX) | (Cnp1 > limI & dXn <-limX) );

%% dXn+1 = f(dXn)
<<<<<<< HEAD

=======
%--- spont ---
[fig, ~] = Spont.a.autocorr_lag_dt(dXissbiais,1)
hold on

%
>>>>>>> master
[binvalsA, elts_per_binA, XnbinReinf] = BinsWithEqualNbofElements(dXn(reinf), dXnp1(reinf), bins1, bins1+3);
[binvalsB, elts_per_binB, XnbinInhib] = BinsWithEqualNbofElements(dXn(inhib), dXnp1(inhib), bins1, bins1+3);
mprodReinf = nanmean(XnbinReinf,2);
mprodInhib = nanmean(XnbinInhib,2);

%***
<<<<<<< HEAD
%get spontaneous trace
[fig, ~] = Spont.autocorr_lag_dt(dXissbiais,1)
hold on
=======
>>>>>>> master
errorbar(binvalsA, mprodReinf, std(XnbinReinf,1,2)/sqrt(elts_per_binA),...
    'Linewidth', 2, 'Color', col_reinf, 'DisplayName', 'reinforcement')
hold on
errorbar(binvalsB, mprodInhib, std(XnbinInhib,1,2)/sqrt(elts_per_binB),...
    'Linewidth', 2, 'Color', col_inhib, 'DisplayName', 'reinforcement')
xlabel('d\theta_n')
ylabel('d\theta_{n+1}')
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 14;

<<<<<<< HEAD
%% C1 = <dXndXn+1>/|dXn||dXn+1|

switch contrast
    case 'relative'
        Contrast_np1 = Cnp1;
    case 'absolute'
        Contrast_np1 = abs(Cnp1);
end

[binvalsC, elts_per_binC, ACbinReinf] = BinsWithEqualNbofElements(Contrast_np1(reinf), ACn_np1(reinf), bins2, bins2);
[binvalsD, elts_per_binD, ACbinInhib] = BinsWithEqualNbofElements(Contrast_np1(inhib), ACn_np1(inhib), bins2, bins2);
=======
%% <dXndXn+1> et <dXndXn+1>/|dXn||dXn+1|

[binvalsA, elts_per_binA, prodbinReinf] = BinsWithEqualNbofElements(Cnp1(reinf), Prodn_np1(reinf), bins2, bins2);
[binvalsB, elts_per_binB, prodbinInhib] = BinsWithEqualNbofElements(Cnp1(inhib), Prodn_np1(inhib), bins2, bins2);
[binvalsC, elts_per_binC, ACbinReinf] = BinsWithEqualNbofElements(Cnp1(reinf), ACn_np1(reinf), bins2, bins2);
[binvalsD, elts_per_binD, ACbinInhib] = BinsWithEqualNbofElements(Cnp1(inhib), ACn_np1(inhib), bins2, bins2);
mprodReinf = nanmean(prodbinReinf,2);
mprodInhib = nanmean(prodbinInhib,2);
>>>>>>> master
mACReinf = nanmean(ACbinReinf,2);
mACInhib = nanmean(ACbinInhib,2);

if forceonzero
<<<<<<< HEAD
=======
    [~, i] = min(abs(binvalsA));
    binvalsA(i) = 0;
    [~, i] = min(abs(binvalsB));
    binvalsB(i) = 0;
>>>>>>> master
    [~, i] = min(abs(binvalsC));
    binvalsC(i) = 0;
    [~, i] = min(abs(binvalsD));
    binvalsD(i) = 0;
end

%***
figure
hold on
<<<<<<< HEAD
errorbar(binvalsC, mACReinf, std(ACbinReinf,1,2)/sqrt(elts_per_binC),...
    '-', 'Linewidth', 2, 'Color', col_reinf, 'DisplayName',...
    '<d\theta_n.d\theta_{n+1}> /<|d\theta_n|.|d\theta_{n+1}|> reinforcement')
errorbar(binvalsD, mACInhib, std(ACbinInhib,1,2)/sqrt(elts_per_binD),...
    '-', 'Linewidth', 2, 'Color', col_inhib, 'DisplayName',...
=======
errorbar(binvalsA, mprodReinf, std(prodbinReinf,1,2)/sqrt(elts_per_binA),...
    'Linewidth', 2, 'Color', col_reinf, 'DisplayName', '<d\theta_n.d\theta_{n+1}> reinforcement')
errorbar(binvalsB, mprodInhib, std(prodbinInhib,1,2)/sqrt(elts_per_binB),...
    'Linewidth', 2, 'Color', col_inhib, 'DisplayName', '<d\theta_n.d\theta_{n+1}> inhibition')
errorbar(binvalsC, mACReinf, std(ACbinReinf,1,2)/sqrt(elts_per_binC),...
    '--', 'Linewidth', 2, 'Color', col_reinf, 'DisplayName',...
    '<d\theta_n.d\theta_{n+1}> /<|d\theta_n|.|d\theta_{n+1}|> reinforcement')
errorbar(binvalsD, mACInhib, std(ACbinInhib,1,2)/sqrt(elts_per_binD),...
    '--', 'Linewidth', 2, 'Color', col_inhib, 'DisplayName',...
>>>>>>> master
    '<d\theta_n.d\theta_{n+1}> /<|d\theta_n|.|d\theta_{n+1}|> inhibition')

legend
xlabel('C_n_+_1 = (I_L-I_R)_n_+_1')
<<<<<<< HEAD
ylabel({'C1'});%, '---', '|(\delta\theta_n*\delta\theta_n_+_1)|'})
=======
ylabel({'<\delta\theta_n*\delta\theta_n_+_1>'});%, '---', '|(\delta\theta_n*\delta\theta_n_+_1)|'})
>>>>>>> master
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 14;


%% Pflip
<<<<<<< HEAD

flips = 0.5*(1 - Prodn_np1./(0.5*Absn_np1));
[binvalsA, elts_per_binA, pflipReinf] = BinsWithEqualNbofElements(Contrast_np1(reinf), flips(reinf), bins2, bins2);
[binvalsB, elts_per_binB, pflipInhib] = BinsWithEqualNbofElements(Contrast_np1(inhib), flips(inhib), bins2, bins2);
=======
[binvalsA, elts_per_binA, absbinReinf] = BinsWithEqualNbofElements((Cnp1(reinf)), Absn_np1(reinf), bins2, bins2);
[binvalsB, elts_per_binB, absBinInhib] = BinsWithEqualNbofElements((Cnp1(inhib)), Absn_np1(inhib), bins2, bins2);
mabsReinf = nanmean(absbinReinf,2);
mabsInhib = nanmean(absBinInhib,2);

p_flip_reinf = 0.5*(1 - mprodReinf./(0.72*mabsReinf));
p_flip_inhib =  0.5*(1 - mprodInhib./(0.72*mabsInhib));

p_flip_reinf_sem = (std(prodbinReinf,1,2)+1./std(absbinReinf,1,2)) /sqrt(elts_per_binA);
p_flip_inhib_sem = (std(prodbinInhib,1,2)+1./std(absBinInhib,1,2)) /sqrt(elts_per_binB);
>>>>>>> master

if forceonzero
    [~, i] = min(abs(binvalsA));
    binvalsA(i) = 0;
    [~, i] = min(abs(binvalsB));
    binvalsB(i) = 0;
end

%***
figure
<<<<<<< HEAD
errorbar(binvalsA, mean(pflipReinf,2), std(pflipReinf,1,2)/sqrt(elts_per_binA),...
    'Linewidth', 2, 'DisplayName', 'reinforcement', 'Color', col_reinf)
hold on
errorbar(binvalsB, mean(pflipInhib,2), std(pflipInhib,1,2)/sqrt(elts_per_binB),...
    'Linewidth', 2, 'DisplayName', 'inhibition', 'Color', col_inhib)
=======
errorbar(binvalsA, p_flip_reinf, p_flip_reinf_sem,...
    'Linewidth', 2, 'DisplayName', 'reinforcement', 'Color', col_reinf)
hold on
errorbar(binvalsB, p_flip_inhib, p_flip_inhib_sem,...
    'Linewidth', 2, 'DisplayName', 'inhibition', 'Color', col_inhib)
%text(-0.6, 0.9, {['pturn = ' num2str(pturn)],['wturn = ' num2str(wturn)]})
>>>>>>> master
legend
ylabel('p_f_l_i_p')
xlabel('C_n_+_1')
ylim([0 1])
ax = gca;
ax.FontName = 'Times New Roman';
<<<<<<< HEAD
ax.FontSize = 20;
=======
ax.FontSize = 14;
>>>>>>> master
