function[] = auto_co_reinforcement(dX, Dii, bins)

colour = colour_palette(0,3);
col_reinf = colour(3,:);
col_inhib = colour(5,:);

Cn = Dii(:, 1:end);
dXn = dX(:,1:end-1);
dXnp1 = dX(:,2:end);
Prodn_np1 = (dXn.*dXnp1);
Absn_np1 = (abs(dXn).*abs(dXnp1));
ACn_np1 = Prodn_np1./Absn_np1;

limI = 0;
limX = 0;
reinf = find( (Cn > limI & dXn >limX) | (Cn < limI & dXn <-limX) );
inhib = find( (Cn < limI & dXn >limX) | (Cn > limI & dXn <-limX) );
%% absolute contrast : <dXndXn+1> et <dXndXn+1>/|dXn||dXn+1|
[binvalsA, elts_per_binA, prodbinReinf] = BinsWithEqualNbofElements(Cn(reinf), Prodn_np1(reinf), round(bins/2), round(bins/2)+1);
[binvalsB, elts_per_binB, prodbinInhib] = BinsWithEqualNbofElements(Cn(inhib), Prodn_np1(inhib), round(bins/2), round(bins/2)+1);
[binvalsC, elts_per_binC, ACbinReinf] = BinsWithEqualNbofElements(Cn(reinf), ACn_np1(reinf), round(bins/2), round(bins/2));
[binvalsD, elts_per_binD, ACbinInhib] = BinsWithEqualNbofElements(Cn(inhib), ACn_np1(inhib), round(bins/2), round(bins/2));
mprodReinf = nanmean(prodbinReinf,2);
mprodInhib = nanmean(prodbinInhib,2);
mACReinf = nanmean(ACbinReinf,2);
mACInhib = nanmean(ACbinInhib,2);

%***
figure
hold on
errorbar(binvalsA, mprodReinf, std(prodbinReinf,1,2)/sqrt(elts_per_binA),...
    'Linewidth', 2, 'Color', col_reinf, 'DisplayName', '<d\theta_n.d\theta_{n+1}> reinforcement')
errorbar(binvalsB, mprodInhib, std(prodbinInhib,1,2)/sqrt(elts_per_binB),...
    'Linewidth', 2, 'Color', col_inhib, 'DisplayName', '<d\theta_n.d\theta_{n+1}> inhibition')
errorbar(binvalsC, mACReinf, std(ACbinReinf,1,2)/sqrt(elts_per_binC),...
    '--', 'Linewidth', 2, 'Color', col_reinf, 'DisplayName',...
    '<d\theta_n.d\theta_{n+1}> /<|d\theta_n|.|d\theta_{n+1}|> reinforcement')
errorbar(binvalsD, mACInhib, std(ACbinInhib,1,2)/sqrt(elts_per_binD),...
    '--', 'Linewidth', 2, 'Color', col_inhib, 'DisplayName',...
    '<d\theta_n.d\theta_{n+1}> /<|d\theta_n|.|d\theta_{n+1}|> inhibition')

legend
xlabel('dI/I')
ylabel({'<\delta\theta_n*\delta\theta_n_+_1>'});%, '---', '|(\delta\theta_n*\delta\theta_n_+_1)|'})
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 14;

%% Pflip
[binvalsA, elts_per_binA, absbinReinf] = BinsWithEqualNbofElements(Cn(reinf), Absn_np1(reinf), round(bins/2), round(bins/2)+1);
[binvalsB, elts_per_binB, absBinInhib] = BinsWithEqualNbofElements(Cn(inhib), Absn_np1(inhib), round(bins/2), round(bins/2)+1);
mabsReinf = nanmean(absbinReinf,2);
mabsInhib = nanmean(absBinInhib,2);

p_flip_reinf = 0.5*(1 - mprodReinf./(0.72*mabsReinf));
p_flip_inhib =  0.5*(1 - mprodInhib./(0.72*mabsInhib));

p_flip_reinf_sem = (std(prodbinReinf,1,2)+1./std(absbinReinf,1,2)) /sqrt(elts_per_binA);
p_flip_inhib_sem = (std(prodbinInhib,1,2)+1./std(absBinInhib,1,2)) /sqrt(elts_per_binB);

%***
figure
errorbar([binvalsA], p_flip_reinf, p_flip_reinf_sem,...
    'Linewidth', 2, 'DisplayName', 'reinforcement', 'Color', col_reinf)
hold on
errorbar([binvalsB],  p_flip_inhib, p_flip_inhib_sem,...
    'Linewidth', 2, 'DisplayName', 'inhibition', 'Color', col_inhib)
%text(-0.6, 0.9, {['pturn = ' num2str(pturn)],['wturn = ' num2str(wturn)]})
legend
ylabel('p_f_l_i_p')
xlabel('dI/I')
ylim([0 1])
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 14;