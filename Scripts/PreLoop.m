%%% Partie de setup avant la boucle d'opti
%%% Chargement des données de cible (si on veut cibler une autre marche que la capture)
%%% Identification des points cibles X

printflag=1;

OPN = PN;
OPol=Pol;
r1 = R_monde_local;
r2 = R_Pelvis_monde_local;
r3 = R_LFem_ref_local;
r4 = R_LTib_ref_local; 
r5 = R_RFem_ref_local; 
r6 = R_RTib_ref_local;

clear Pol;

% Chargement du polynôme de la spline de la trajectoire angulaire et des paramètres physiologiques à cibler
% On prend le polynôme et non directement la courbe, car le nombre de frame
% par simulation - et donc l'échantillonage - change à chaque fichier...
load(strcat(PathPreSet,Names{jj},'.mat')); 

% Reconstruction de la courbe de TA
NewCurve=[];
% for i =1:Pol(end,1)
%     NewCurve = [NewCurve , spline_to_curve_int(Pol(Pol(:,1)==i,:),Period)'];
% end
% if NewCurve(end,:) ~= NewCurve(1,:)
%     NewCurve(end,:) = NewCurve(1,:);
% end


NewCurve=[];
for i =1:Period+1
    tmp=zeros(1,11);
    for j=1:11
        tmp(j) = EvalSpline(Pol(Pol(:,1)==j,:),((i-1))/Period);
    end
    NewCurve = [NewCurve ; tmp];
end

% NewCurve=OCurve;

% 
PN=[];
for i =1:Period+1
    PN = [PN; fcine_numerique_H2(NewCurve(i,:),Param(:,1)',Param(:,2)',Param(:,3)',Param(:,4)',Param(:,5)',Param(:,6)', R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local)'];
end

% if PN(end,:) ~= PN(1,:)
%     PN(end,:) = PN(1,:);
% end
% Courbe cible devient la nouvelle poulaine
GT = PN;
% GT(:,4:6) = [GT(mid:end-1,1) , -GT(mid:end-1,2) , GT(mid:end-1,3) ; GT(1:mid,1) , -GT(1:mid,2) ,GT(1:mid,3) ];
% PN = GT;

% Rétablissement des anciens paramètres 
% Pol = OPol;
PN = OPN;
Param = [Fem1g',Fem6g',Tal1g', Fem1d', Fem6d', Tal1d'];

R_monde_local = r1;
R_Pelvis_monde_local = r2;
R_LFem_ref_local = r3;
R_LTib_ref_local = r4; 
R_RFem_ref_local = r5; 
R_RTib_ref_local = r6;


%%%%%%%%%%

% Target X Setup :

[~,ag] = min(GT(:,1));
[~,bg] = max(GT(:,1));
[~,cg] = min(GT(:,3));

[~,ad] = min(GT(:,4));
[~,bd] = max(GT(:,4));
[~,cd] = min(GT(:,6));

% A = [ag;bg;cg;mod(ag+mid,Period);mod(bg+mid,Period);mod(cg+mid,Period)];
% A = [ag;bg;cg;ad;bd;cd];
A = [ag;bg;cg];
A = [A ; mod(A+mid+1,Period)+1];
if ~isempty(A==0)
    % 1-périodicité, pour éviter un indice 0 dans une ancienne méthode
    A(A==0) = Period;
end

% Création du vecteur des points cibles - X 
% La première colonne repère quel est la courbe dont est issu le point originalement. C'était pour le calcul des 
% tangentes, mais c'est pas utile actuellement.
X = [2;2;3;5;5;6];
X = [X , (A-1)/Period];
X = [X , [GT(A(1:3),1:3) ; GT(A(4:6),4:6)]];

% X(:,4) = zeros(size(6,1));

% Introduction d'un delta -> tests
% X(1,3) = X(1,3) + 0.1;
% X(4,3) = X(4,3) + 0.1;

% X(1,2) = X(1,2) + 0.1;
% X(4,2) = X(4,2) + 0.1;

% X(1,3) = X(1,3) + 0.05;
% X(4,3) = X(4,3) + 0.05;

% X(1,3) = X(1,3) + 0.01;
% X(4,3) = X(4,3) + 0.01;

if ii==jj
    R = rand(6,3);
    X(:,3:5) = X(:,3:5) + 0.01*(R/norm(R));
end
% Affichage Courbes poulaines Cible GT, Construite à partir des TA PN, Issue de la capture splinée NewP
% et les cibles X
if printflag
figure;
hold on;
for i =1:6
    subplot(2,3,i);
    hold on;
%     plot(1:size(P,1),P(:,i),'b');
%     plot(1:size(PN,1),PN(:,i),'g');
    plot(0:1/Period:1,PN(:,i),'g');
%     plot(1:size(NewP), NewP(:,i),'r');
%     plot(1:size(GT), GT(:,i),'k');
    plot(0:1/Period:1, GT(:,i),'k');
    if i <4
%         plot(X(1:3,2)*Period, X(1:3,i+2),'kx');
        plot(X(1:3,2), X(1:3,i+2),'kx');
    else
%         plot(X(4:6,2)*Period, X(4:6,i-1),'kx');
        plot(X(4:6,2), X(4:6,i-1),'kx');
    end
end
norm(X(:,3:5)-ErrorFun(OPol,X,Param, R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local))
norm(X(:,3:5)-ErrorFun2(Pol,X,Param, R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local))
norm(X(:,3:5)-ErrorFun2(OPol,X,Param, R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local))

end

PN;
