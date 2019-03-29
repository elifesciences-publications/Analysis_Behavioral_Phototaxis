function f = double_gauss_constrained(x, a, alpha, beta)

%alpha=0.2701;
%beta=0.1727;

w1 = (alpha*a+sqrt(alpha^2*a^2-(alpha^2-beta*(1-a))*(a^2+a*(1-a))))/(a^2+a*(1-a));
w2 = (alpha*sqrt(pi/2)-a*w1)/(1-a);
a2 = 1-a;
disp([a, a2, w1, w2])

f = a/sqrt(2*pi)/w1*exp(-(x/w1).^2/2)+a2/sqrt(2*pi)/w2*exp(-(x/w2).^2/2);


end

