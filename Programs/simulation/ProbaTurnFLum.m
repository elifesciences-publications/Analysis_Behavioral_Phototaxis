function[p_turn, wturn] = ProbaTurnFLum(dll)

basal_wturn = 0.59; 
basal_pturn = 0.41;

p_turn = -dll*0.5 + basal_pturn;
p_turn(dll>0) = basal_pturn;
p_turn(dll>0.3) = basal_pturn+0.1;
p_turn(dll>0.5 & dll<1) = basal_pturn+0.2;
p_turn(p_turn>1) = 1;

wturn = -2*dll + basal_wturn;
%wturn = -2*dll + basal_wturn;
wturn(dll>0) = basal_wturn-0.3*dll(dll>0);
 
end