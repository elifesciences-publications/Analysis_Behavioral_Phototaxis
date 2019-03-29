function[p_turn, wturn] = ProbaTurnFLum(dll)

basal_wturn = 0.6; 
basal_pturn = 0.5;

p_turn = -dll*0.33 + basal_pturn;
p_turn(dll>0) = basal_pturn;
p_turn(abs(dll)<0.2) = 0.4;

wturn = -dll + basal_wturn;
wturn(dll>0) = basal_wturn;
 
end