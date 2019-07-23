function[f] = fit_autocorrIBI(x, rflip)

N=10000000;
trans2 = [1-rflip,rflip;rflip,1-rflip];
emis2=[0,1; 1,0];
[seq2,~] = hmmgenerate(N,trans2,emis2);

 
% Numerical generation of scoot and turn states
pturn=0.41;        % fraction of turning events
pt2f=0.47;         % this value is set by dx2n-1 vs dx2n fitting
pf2t=pt2f*pturn/(1-pturn);
pt2f=1-pturn;
pf2t=pturn;
trans3=[1-pt2f,pt2f;pf2t,1-pf2t];
emis3=[0,1; 1,0];
[seq3,~] = hmmgenerate(N,trans3,emis3);


% generate angular sequence
wturn=0.59;  % std of the turn distribution in radian
wfor=0.092;   % sdt of the forward distribution in radian

 
seq=[seq1;seq2;seq3];
Tbouts=find(seq1==2);
%seq_bouts=seq(:,Tbouts);
ind_turn=find(seq(3,:)==2 & seq(1,:)==2);
ind_for=find(seq(3,:)==1 & seq(1,:)==2);
angle_turns=randn(1,length(ind_turn))*wturn;
angle_for=randn(1,length(ind_for))*wfor;

 
angles=zeros(N,1);
angles(ind_for)=angle_for;
sign_angle_turn=[seq2(ind_turn)>1.5]*2-1;
angles(ind_turn)=sign_angle_turn.*abs(angle_turns);
angles_disc=angles(Tbouts);

a=autocorr(angles_disc, 20);


end