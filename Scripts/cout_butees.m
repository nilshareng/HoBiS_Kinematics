function res=cout_butees(angle,bmin,bmax)
% fonction en exponentiel sur les but�es

coef=100;
% gestion de la but�e max
res=exp(coef*(angle-bmax));

% gestion de la but�e min
res=res+exp(coef*(bmin-angle));
