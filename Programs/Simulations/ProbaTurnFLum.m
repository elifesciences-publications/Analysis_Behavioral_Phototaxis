function[p_turn, wturn] = ProbaTurnFLum(dii)

basal_wturn = 0.59; 
basal_pturn = 0.41;

%%
% --- pturn ---
% negative dii
p_turn = -dii*0.1 + 0.0513 + basal_pturn;
%plateau
%p_turn(dii>-0.12) = 0.3*dii(dii>-0.12) + 0.05 + basal_pturn;
p_turn(dii>-0.2 & dii<0.2) = basal_pturn;
p_turn(dii>0.2) = 0.5;
% upper threshold
p_turn(p_turn>1) = 1;

% --- wturn ---
wturn = -0.5*dii + basal_wturn; %0.72
wturn(dii>0) = -0.1*dii(dii>0) + basal_wturn;

%%
p_turn = ones(size(dii))*basal_pturn;
wturn = ones(size(dii))*basal_wturn;
 
end