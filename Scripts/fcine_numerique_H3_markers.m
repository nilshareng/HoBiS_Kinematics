function [res,Save,Nmarkers]=fcine_numerique_H3_markers(P,markers,Fem1g,Fem6g,Tal1g,Fem1d,Fem6d,Tal1d, R_monde_local, R_Pelvis_monde, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local)
% Fonction cin�matique du mod�le
% Entr�e : Angles P et param�tres physiologiques Fem...
% Param en ligne
% Sortie : Les positions des chevilles simul�es
% Franck : int�gration des rotations anatomiques initiales

newmarkers = markers;

% Parseur
rxb=P(1);
ryb=P(2);
rzb=P(3);
rhxg=P(4);
rhyg=P(5);
rhzg=P(6); 
rzgg=P(7);
rhxd=P(8);
rhyd=P(9);
rhzd=P(10); 
rzgd=P(11);

flag=0;
if size(P,2) == 12
    flag =1;
end

Mpass = eye(4); % Monde 

Mpass = Mpass * R_Pelvis_monde; % On se pose en Opelvis, orientation pelvis de reference   

% Rotation angles pelvis
matrbasX=[1 0 0 0;  0 cos(rxb) -sin(rxb) 0; 0 sin(rxb) cos(rxb) 0 ; 0 0 0 1];
matrbasY=[cos(ryb) 0 sin(ryb) 0;  0 1 0 0; -sin(ryb) 0 cos(ryb) 0 ; 0 0 0 1];
matrbasZ=[cos(rzb) -sin(rzb) 0 0;  sin(rzb) cos(rzb) 0 0; 0 0 1 0; 0 0 0 1];
matrbas=matrbasX*matrbasY*matrbasZ;

% multiplication par R_pelvis_monde pour passer du rep�re du monde (rappel
% : Y_up (0;1;0) parfait, X vers l'avant, Z lat�ral) vers le rep�re du
% bassin l�g�rement inclin� par rapport aux axes du monde
 
% temp=eye(4,4); temp(1:3,1:3)=R_Pelvis_monde;
Mpass= Mpass * matrbas; 

% temp(1:3,1:3)=R_monde_local; % ramener au rep�re du monde avec Z en haut
% Mpass= temp*Mpass; 

%actualisation de matrbas qui va �tre r�utilis� pour le cot� droit
matrbas=Mpass;

%exprimer Fem1g dans le rep�re du bassin, alors qu'il est initialement
%exprim� dans le rep�re du monde positionn� sur le bassin : R_Pelvis_monde
%(pelvis vers rep�re local monde et repere local_monde vers monde avec Z-up)
% Fem1glocal=(R_monde_local*R_Pelvis_monde)'*Fem1g'; 
% Fem1glocal = Mpass(1:3,1:3) * Fem1g;

% mattbth=[1 0 0 Fem1glocal(1) ; 0 1 0 Fem1glocal(2) ;0 0 1 Fem1glocal(3) ; 0 0 0 1];
Mpass= Mpass*R_LFem_ref_local;
Save = [];
Save = [Save , 1;Mpass(1:3,4)];


%Franck : insertion de la matrice pour aligner sur posture de
%ref/description
% mathancheg=eye(4,4);
% mathancheg(1:3,1:3)=R_LFem_ref_local; 

%matrice de rotation li�e au mouvement hanche gauche
mathancheXg=[1 0 0 0;  0 cos(rhxg) -sin(rhxg) 0; 0 sin(rhxg) cos(rhxg) 0 ; 0 0 0 1];
mathancheYg=[cos(rhyg) 0 sin(rhyg) 0;  0 1 0 0; -sin(rhyg) 0 cos(rhyg) 0 ; 0 0 0 1];
mathancheZg=[cos(rhzg) -sin(rhzg) 0 0;  sin(rhzg) cos(rhzg) 0 0; 0 0 1 0; 0 0 0 1];

if flag
    mathancheg=mathancheZg*mathancheYg*mathancheXg;
else
    mathancheg=mathancheZg*mathancheYg*mathancheXg;
    %mathancheg=mathancheXg*mathancheYg*mathancheZg;
end

Mpass=Mpass*mathancheg; 

% Zp=Fem6g-Fem1g; 
% Zp=(R_monde_local*R_Pelvis_monde*R_LFem_ref_local)'*Zp'; 


% matF1F6g=[1 0 0 Zp(1) ; 0 1 0 Zp(2) ;0 0 1 Zp(3) ; 0 0 0 1];
% matF1F6g=[1 0 0 0 ; 0 1 0 -norm(Zp) ;0 0 1 0 ; 0 0 0 1];

Mpass=Mpass*R_LTib_ref_local;

Save = [Save , [2;Mpass(1:3,4)]];

%Franck : alignement du tibia avec la posture de ref/description

%matzgg=[1 0 0 0;  0 cos(rzgg) -sin(rzgg) 0; 0 sin(rzgg) cos(rzgg) 0 ; 0 0 0 1];
matzgg=[cos(rzgg) -sin(rzgg) 0 0; sin(rzgg) cos(rzgg) 0 0 ; 0 0 1 0; 0 0 0 1];

if flag
    matzgg=[cos(rzgg) -sin(rzgg) 0 0;  sin(rzgg) cos(rzgg) 0 0; 0 0 1 0; 0 0 0 1];
end
% temp=eye(4,4); temp(1:3,1:3)=R_LTib_ref_local;
Mpass=Mpass*matzgg;

ttz=Tal1g;
MtransiTib = (R_Pelvis_monde * R_LFem_ref_local * R_LTib_ref_local)^-1;
newmarkers.RKNI
newmarkers.LKNI
newmarkers.RANE
newmarkers.RANI
newmarkers.LANE
newmarkers.LANI

newmarkers.RANE
newmarkers.RANI
newmarkers.LANE
newmarkers.LANI
;
ttz = MtransiTib * [ttz';1];
% ttz=(R_Pelvis_monde*R_LFem_ref_local^-1*R_LTib_ref_local)'*ttz'; % A modif car Fem & Tib non ortho dir 
% MPass = MPass(1:3,1:3) *  ttz';
% matF6T1g=[1 0 0 ttz(1) ; 0 1 0 ttz(2) ;0 0 1 ttz(3) ; 0 0 0 1];
% matF6T1g=[1 0 0 0 ; 0 1 0 -norm(ttz) ;0 0 1 0 ; 0 0 0 1];
% Mpass=Mpass*matF6T1g;
% Save = [Save , [3;Mpass(1:3,4)]];
ttz = Mpass *  ttz;
Save = [Save , [3;ttz(1:3)]];




% Bilan : Fonction de cinematique directe : 
matg=ttz(1:3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% DROITE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mpass= matrbas;

%exprimer Fem1d dans le rep�re du bassin, alors qu'il est initialement
%exprim� dans le rep�re du monde positionn� sur le bassin : R_Pelvis_monde
%(pelvis vers rep�re local monde) puis passage au rep�re du monde Z-UP
%Fem1dlocal=R_Pelvis_monde'*Fem1d'; 
% Fem1dlocal=(R_monde_local*R_Pelvis_monde)'*Fem1d'; 



% mattbth=[1 0 0 Fem1dlocal(1) ; 0 1 0 Fem1dlocal(2) ;0 0 1 Fem1dlocal(3) ; 0 0 0 1];
Mpass= Mpass*R_RFem_ref_local;


Save = [Save , [4;Mpass(1:3,4)]];

%Franck : insertion de la matrice pour aligner sur posture de
%ref/description
% mathanched=eye(4,4);
% mathanched(1:3,1:3)=R_RFem_ref_local; 

mathancheXd=[1 0 0 0;  0 cos(rhxd) -sin(rhxd) 0; 0 sin(rhxd) cos(rhxd) 0 ; 0 0 0 1];
mathancheYd=[cos(rhyd) 0 sin(rhyd) 0;  0 1 0 0; -sin(rhyd) 0 cos(rhyd) 0 ; 0 0 0 1];
mathancheZd=[cos(rhzd) -sin(rhzd) 0 0;  sin(rhzd) cos(rhzd) 0 0; 0 0 1 0; 0 0 0 1];

if flag
    mathanched=mathancheZd*mathancheYd*mathancheXd;
else
    %mathanched=mathancheXd*mathancheYd*mathancheZd;
    mathanched=mathancheZd*mathancheYd*mathancheXd;
end

Mpass=Mpass*mathanched; 

% Zp=Fem6d-Fem1d; 
% Zp=(R_monde_local*R_Pelvis_monde*R_RFem_ref_local)'*Zp'; 

% matF1F6d=[1 0 0 Zp(1) ; 0 1 0 Zp(2) ;0 0 1 Zp(3) ; 0 0 0 1];
% matF1F6d=[1 0 0 0 ; 0 1 0 -norm(Zp) ;0 0 1 0 ; 0 0 0 1];

% Mpass=Mpass*matF1F6d;
Mpass = Mpass * R_RTib_ref_local;

Save = [Save , [5;Mpass(1:3,4)]];


%rzfemg=0;
%mattorsfemg=[cos(rzfemg) -sin(rzfemg) 0 0;  sin(rzfemg) cos(rzfemg) 0 0; 0 0 1 0; 0 0 0 1];

%Mpass=Mpass*mattorsfemg;

% Matrice de passage pour rendre l'axe du genou (fem9-fem10) parrallele � l'axe du tibia (tib5-tib6)
%rygeng=0;
%Mpass_genoug=[cos(rygeng) 0 sin(rygeng) 0;  0 1 0 0; -sin(rygeng) 0 cos(rygeng) 0 ; 0 0 0 1];

%Mpass=Mpass*Mpass_genoug;
%matzgd=[1 0 0 0;  0 cos(rzgd) -sin(rzgd) 0; 0 sin(rzgd) cos(rzgd) 0 ; 0 0 0 1];
matzgd=[cos(rzgd) -sin(rzgd) 0 0; sin(rzgd) cos(rzgd) 0 0; 0 0 1 0; 0 0 0 1];

if flag
    matzgd=[cos(rzgd) -sin(rzgd) 0 0;  sin(rzgd) cos(rzgd) 0 0; 0 0 1 0; 0 0 0 1];
end

% temp=eye(4,4); temp(1:3,1:3)=R_RTib_ref_local;
% Mpass=Mpass*temp*matzgd;
Mpass=Mpass*matzgd;

ttz = Tal1d;
ttz = (R_Pelvis_monde * R_RFem_ref_local * R_RTib_ref_local)^-1 * [ttz';1];
ttz = Mpass *  ttz;

Save = [Save , [6;ttz(1:3)]];
% ttz=(R_monde_local*R_Pelvis_monde*R_RFem_ref_local*R_RTib_ref_local)'*ttz'; 

% matF6T1d=[1 0 0 ttz(1) ; 0 1 0 ttz(2) ;0 0 1 ttz(3) ; 0 0 0 1];
% matF6T1d=[1 0 0 0 ; 0 1 0 -norm(ttz) ;0 0 1 0 ; 0 0 0 1];
% Mpass=Mpass*matF6T1d;

% Save = [Save , [6;ttz(1:3)]];

% Bilan : Fonction de cinematique directe : 
% matd=Mpass*[0;0;0;1];
matd=ttz(1:3);
% R�organisation en fonction du sens de la marche, rep�r� par la coordonn�e lat�rale de la Hanche 
% if (Fem1g(1)>0)
%     res = [matd(1:2) ; matg(3) ; matg(1:2) ; matd(3)];
% else
    res= [matg(1:3);matd(1:3)];
% end

