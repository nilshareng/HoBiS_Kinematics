function res=cout_butees(angle,bmin,bmax)
% fonction en exponentiel sur les butées

coef=100;
% gestion de la butée max
res=exp(coef*(angle-bmax));

% gestion de la butée min
res=res+exp(coef*(bmin-angle));
