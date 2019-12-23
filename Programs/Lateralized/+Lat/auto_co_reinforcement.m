function[] = auto_co_reinforcement(dXl,dXissbiais, DIlr, bins1, bins2, forceonzero, contrast)

% --- appearance ---
colour = colour_palette(0,3);
col_reinf = colour(3,:);
col_inhib = colour(4,:);

% --- 
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

[binvalsA, elts_per_binA, XnbinReinf] = BinsWithEqualNbofElements(dXn(reinf), dXnp1(reinf), bins1, bins1+3);
[binvalsB, elts_per_binB, XnbinInhib] = BinsWithEqualNbofElements(dXn(inhib), dXnp1(inhib), bins1, bins1+3);
mprodReinf = nanmean(XnbinReinf,2);
mprodInhib = nanmean(XnbinInhib,2);

%***
%get spontaneous trace
[fig, ~] = Spont.autocorr_lag_dt(dXissbiais,1)
hold on
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

%% C1 = <dXndXn+1>/|dXn||dXn+1|

switch contrast
    case 'relative'
        Contrast_np1 = Cnp1;
    case 'absolute'
        Contrast_np1 = abs(Cnp1);
end

[binvalsC, elts_per_binC, ACbinReinf] = BinsWithEqualNbofElements(Contrast_np1(reinf), ACn_np1(reinf), bins2, bins2);
[binvalsD, elts_per_binD, ACbinInhib] = BinsWithEqualNbofElements(Contrast_np1(inhib), ACn_np1(inhib), bins2, bins2);
mACReinf = nanmean(ACbinReinf,2);
mACInhib = nanmean(ACbinInhib,2);

if forceonzero
    [~, i] = min(abs(binvalsC));
    binvalsC(i) = 0;
    [~, i] = min(abs(binvalsD));
    binvalsD(i) = 0;
end

%***
figure
hold on
errorbar(binvalsC, mACReinf, std(ACbinReinf,1,2)/sqrt(elts_per_binC),...
    '-', 'Linewidth', 2, 'Color', col_reinf, 'DisplayName',...
    '<d\theta_n.d\theta_{n+1}> /<|d\theta_n|.|d\theta_{n+1}|> reinforcement')
errorbar(binvalsD, mACInhib, std(ACbinInhib,1,2)/sqrt(elts_per_binD),...
    '-', 'Linewidth', 2, 'Color', col_inhib, 'DisplayName',...
    '<d\theta_n.d\theta_{n+1}> /<|d\theta_n|.|d\theta_{n+1}|> inhibition')

legend
xlabel('C_n_+_1 = (I_L-I_R)_n_+_1')
ylabel({'C1'});%, '---', '|(\delta\theta_n*\delta\theta_n_+_1)|'})
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 14;


%% Pflip

%flips = 0.5*(1 - Prodn_np1./(0.65*Absn_np1));
flips = 0.5*(1 - 1/0.41*Prodn_np1./Absn_np1);

[binvalsA, elts_per_binA, pflipReinf] = BinsWithEqualNbofElements(Contrast_np1(reinf), flips(reinf), bins2, bins2);
[binvalsB, elts_per_binB, pflipInhib] = BinsWithEqualNbofElements(Contrast_np1(inhib), flips(inhib), bins2, bins2);

if forceonzero
    [~, i] = min(abs(binvalsA));
    binvalsA(i) = 0;
    [~, i] = min(abs(binvalsB));
    binvalsB(i) = 0;
end

%***
figure
errorbar(binvalsA, mean(pflipReinf,2), std(pflipReinf,1,2)/sqrt(elts_per_binA),...
    'Linewidth', 2, 'DisplayName', 'reinforcement', 'Color', col_reinf)
hold on
errorbar(binvalsB, mean(pflipInhib,2), std(pflipInhib,1,2)/sqrt(elts_per_binB),...
    'Linewidth', 2, 'DisplayName', 'inhibition', 'Color', col_inhib)
legend
ylabel('p_f_l_i_p')
xlabel('C_n_+_1')
ylim([0 1])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 20;
