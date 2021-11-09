function er = energie_cinetique(etatc,etatt,M, Fem1g, Fem1d, Fem6g, Fem6d, Tal1g, Tal1d,poulainegc, poulainegt,poulainedc, poulainedt, appui,Inertie)
% fonction qui calcule le w en fonction de l'état courant (avant IK) et de l'état à t(après IK)
% en utilisant la matrice de préproduit vectoriel

dt=1/60;

matrbasX=[1 0 0;  0 cos(etatc(1)) -sin(etatc(1)); 0 sin(etatc(1)) cos(etatc(1)) ];
matrbasY=[cos(etatc(2)) 0 sin(etatc(2)) ;  0 1 0 ; -sin(etatc(2)) 0 cos(etatc(2)) ];
matrbasZ=[cos(etatc(3)) -sin(etatc(3)) 0 ;  sin(etatc(3)) cos(etatc(3)) 0 ; 0 0 1 ];
matrbasc=matrbasX*matrbasY*matrbasZ;
matrbasX=[1 0 0;  0 cos(etatt(1)) -sin(etatt(1)); 0 sin(etatt(1)) cos(etatt(1)) ];
matrbasY=[cos(etatt(2)) 0 sin(etatt(2)) ;  0 1 0 ; -sin(etatt(2)) 0 cos(etatt(2)) ];
matrbasZ=[cos(etatt(3)) -sin(etatt(3)) 0 ;  sin(etatt(3)) cos(etatt(3)) 0 ; 0 0 1 ];
matrbast=matrbasX*matrbasY*matrbasZ;

% calcul de l'inertie du bassin+tronc
R=matrbasc;
dR = (matrbast - matrbasc)/(dt);
Mw = dR*R';
w(1) = Mw(3,2);
w(2) = Mw(1,3);
w(3) = Mw(2,1);
if nargin == 15
    I = Inertie.Pelvis;
else
    I = inertie(1.70,50*M,1);
end% pour une taille de 1.70m , une masse de 60Kg et le segment tronc (1)
er = norm(0.5.*((R*I*R'*(w(:).^2)))'); 

Mpassc= matrbasc;
Mpasst= matrbast;

pFem1gc=[[Mpassc zeros(3,1)]; [0 0 0 1]]*[1 0 0 Fem1g(1) ; 0 1 0 Fem1g(2) ;0 0 1 Fem1g(3) ; 0 0 0 1];
pFem1gt=[[Mpasst zeros(3,1)]; [0 0 0 1]]*[1 0 0 Fem1g(1) ; 0 1 0 Fem1g(2) ;0 0 1 Fem1g(3) ; 0 0 0 1];

mathancheXg=[1 0 0 0;  0 cos(etatc(4)) -sin(etatc(4)) 0; 0 sin(etatc(4)) cos(etatc(4)) 0; 0 0 0 1];
mathancheYg=[cos(etatc(5)) 0 sin(etatc(5)) 0;  0 1 0 0; -sin(etatc(5)) 0 cos(etatc(5)) 0; 0 0 0 1];
mathancheZg=[cos(etatc(6)) -sin(etatc(6)) 0 0;  sin(etatc(6)) cos(etatc(6)) 0 0; 0 0 1 0; 0 0 0 1];
mathanchegc=mathancheXg*mathancheYg*mathancheZg;
mathancheXg=[1 0 0 0;  0 cos(etatt(4)) -sin(etatt(4)) 0; 0 sin(etatt(4)) cos(etatt(4)) 0; 0 0 0 1];
mathancheYg=[cos(etatt(5)) 0 sin(etatt(5)) 0;  0 1 0 0; -sin(etatt(5)) 0 cos(etatt(5)) 0; 0 0 0 1];
mathancheZg=[cos(etatt(6)) -sin(etatt(6)) 0 0;  sin(etatt(6)) cos(etatt(6)) 0 0; 0 0 1 0; 0 0 0 1];
mathanchegt=mathancheXg*mathancheYg*mathancheZg;

Mpassc=pFem1gc*mathanchegc; 
Mpasst=pFem1gt*mathanchegt; 
Rc=Mpassc(1:3,1:3);
dR = (Mpasst(1:3,1:3) - Mpassc(1:3,1:3))/(dt);
Mw = dR*R';
wc(1) = Mw(3,2);
wc(2) = Mw(1,3);
wc(3) = Mw(2,1);
if nargin == 15
    Ic = Inertie.Hip;
else
    Ic= inertie(1.70,M,2); % pour une taille de 1.70m , une masse de 60Kg et le segment cuisse (1)
end
% il faut rajouter l'inertie du tibia
%er = er+norm(0.5.*((R*I*R'*(w(:).^2)))'); 

Zp=(Fem6g-Fem1g); 
cuissegc=Mpassc*[0.41*Zp(1) ; 0.41*Zp(2) ;0.41*Zp(3) ;  1];
cuissegt=Mpasst*[0.41*Zp(1) ; 0.41*Zp(2) ;0.41*Zp(3) ;  1];
dcuisseg=(cuissegt(1:3)-cuissegc(1:3))/dt;

matF1F6g=[1 0 0 Zp(1) ; 0 1 0 Zp(2) ;0 0 1 Zp(3) ; 0 0 0 1];    

Mpassc=Mpassc*matF1F6g; 
Mpasst=Mpasst*matF1F6g; 
genout=Mpasst*[0;0;0;1];

matxggc=[1 0 0 0;  0 cos(etatc(7)) -sin(etatc(7)) 0; 0 sin(etatc(7)) cos(etatc(7)) 0; 0 0 0 1];
matxggt=[1 0 0 0;  0 cos(etatt(7)) -sin(etatt(7)) 0; 0 sin(etatt(7)) cos(etatt(7)) 0; 0 0 0 1];
Mpassc=Mpassc*matxggc;
Mpasst=Mpasst*matxggt;
Zp=(Tal1g-Fem6g); 
tibiagc=Mpassc*[Zp(1);Zp(2);Zp(3);1];
tibiagt=Mpasst*[Zp(1);Zp(2);Zp(3);1];   
dtibiag=(tibiagt(1:3)-tibiagc(1:3))/dt;

Rt=Mpassc(1:3,1:3);
dR = (Mpasst(1:3,1:3) - Mpassc(1:3,1:3))/(dt);
Mw = dR*R';
wt(1) = Mw(3,2);
wt(2) = Mw(1,3);
wt(3) = Mw(2,1);
if nargin == 15
    It = Inertie.Leg;
else
    It = inertie(1.70,M,3); % pour une taille de 1.70m , une masse de 60Kg et le segment tibia 3)
end
temp=0.433*norm((genout(1:3)-pFem1gt(1:3,4)));
pmasse = [0.497,0.1,0.0465];
Itemp=eye(3,3); Itemp(1,1)=temp^2; Itemp(2,2)=temp^2; Itemp(3,3)=temp^2;
Itemp=Ic+It+masse(M,3,pmasse)*Itemp;
er = er+norm(0.5.*((Rc*Itemp*Rc'*(wc(:).^2)))'); % cuisse + tibia
er = er+norm(0.5.*((Rt*It*Rt'*(wt(:).^2)))'); % tibia isolément


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% DROITE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mpassc= matrbasc;
Mpasst= matrbast;

pFem1gc=[[Mpassc zeros(3,1)]; [0 0 0 1]]*[1 0 0 Fem1d(1) ; 0 1 0 Fem1d(2) ;0 0 1 Fem1d(3) ; 0 0 0 1];
pFem1gt=[[Mpasst zeros(3,1)]; [0 0 0 1]]*[1 0 0 Fem1d(1) ; 0 1 0 Fem1d(2) ;0 0 1 Fem1d(3) ; 0 0 0 1];

mathancheXg=[1 0 0 0;  0 cos(etatc(8)) -sin(etatc(8)) 0; 0 sin(etatc(8)) cos(etatc(8)) 0; 0 0 0 1];
mathancheYg=[cos(etatc(9)) 0 sin(etatc(9)) 0;  0 1 0 0; -sin(etatc(9)) 0 cos(etatc(9)) 0; 0 0 0 1];
mathancheZg=[cos(etatc(10)) -sin(etatc(10)) 0 0;  sin(etatc(10)) cos(etatc(10)) 0 0; 0 0 1 0; 0 0 0 1];
mathanchegc=mathancheXg*mathancheYg*mathancheZg;
mathancheXg=[1 0 0 0;  0 cos(etatt(8)) -sin(etatt(8)) 0; 0 sin(etatt(8)) cos(etatt(8)) 0; 0 0 0 1];
mathancheYg=[cos(etatt(9)) 0 sin(etatt(9)) 0;  0 1 0 0; -sin(etatt(9)) 0 cos(etatt(9)) 0; 0 0 0 1];
mathancheZg=[cos(etatt(10)) -sin(etatt(10)) 0 0;  sin(etatt(10)) cos(etatt(10)) 0 0; 0 0 1 0; 0 0 0 1];
mathanchegt=mathancheXg*mathancheYg*mathancheZg;

Mpassc=pFem1gc*mathanchegc; 
Mpasst=pFem1gt*mathanchegt; 
R=Mpassc(1:3,1:3);
dR = (Mpasst(1:3,1:3) - Mpassc(1:3,1:3))/(dt);
Mw = dR*R';
w(1) = Mw(3,2);
w(2) = Mw(1,3);
w(3) = Mw(2,1);
if nargin == 15
    I = Inertie.Hip;
else
    I = inertie(1.70,M,2); % pour une taille de 1.70m , une masse de 60Kg et le segment tronc (1)
end

er = er+ norm(0.5.*((R*I*R'*(w(:).^2)))'); 

Zp=(Fem6d-Fem1d); 
cuissegc=Mpassc*[0.41*Zp(1) ; 0.41*Zp(2) ;0.41*Zp(3) ;  1];
cuissegt=Mpasst*[0.41*Zp(1) ; 0.41*Zp(2) ;0.41*Zp(3) ;  1];
dcuissed=(cuissegt(1:3)-cuissegc(1:3))/dt;

matF1F6g=[1 0 0 Zp(1) ; 0 1 0 Zp(2) ;0 0 1 Zp(3) ; 0 0 0 1];

Mpassc=Mpassc*matF1F6g; 
Mpasst=Mpasst*matF1F6g; 

matxggc=[1 0 0 0;  0 cos(etatc(11)) -sin(etatc(11)) 0; 0 sin(etatc(11)) cos(etatc(11)) 0; 0 0 0 1];
matxggt=[1 0 0 0;  0 cos(etatt(11)) -sin(etatt(11)) 0; 0 sin(etatt(11)) cos(etatt(11)) 0; 0 0 0 1];
Mpassc=Mpassc*matxggc;
Mpasst=Mpasst*matxggt;
Zp=(Tal1d-Fem6d); 
tibiagc=Mpassc*[Zp(1);Zp(2);Zp(3);1];
tibiagt=Mpasst*[Zp(1);Zp(2);Zp(3);1];
dtibiad=(tibiagt(1:3)-tibiagc(1:3))/dt;

R=Mpassc(1:3,1:3);
dR = (Mpasst(1:3,1:3) - Mpassc(1:3,1:3))/(dt);
Mw = dR*R';
w(1) = Mw(3,2);
w(2) = Mw(1,3);
w(3) = Mw(2,1);
if nargin == 15
    I = Inertie.Leg;
else
    I = inertie(1.70,M,3); % pour une taille de 1.70m , une masse de 60Kg et le segment tronc (1)
end
er = er+norm(0.5.*((R*I*R'*(w(:).^2)))'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mt=M*0.497;
mc=M*0.10;
mtib=M*0.0465;
if (appui==0)
  ind=1:3;
  VG=norm((poulainegt-poulainegc)/dt);
else 
%  if (appui==1)
    ind=4:6;
    VG=norm((poulainedt-poulainedc)/dt);
%   else 
%     ind=1:3;
%     VGg=norm((poulainet(ind)-poulainec(ind))/dt);
%     ind=4:6;
%     VGd=norm((poulainet(ind)-poulainec(ind))/dt);
%     VG=appui*VGd+(1-appui)*VGg;
%     VGt=((1-appui)*(poulainec(1:3)-poulainet(1:3))+(appui)*(poulainec(4:6)-poulainet(4:6)))/dt;
%   end;
end
ec=0.5*mc*(norm(dcuisseg)^2+norm(dcuissed)^2)+0.5*mtib*(norm(dtibiag)^2+norm(dtibiad)^2)+0.5*M*(VG)^2;

er=er+ec;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FONCTION INERTIE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ii = inertie(hauteur,mglobal,membres)

% membres = 1 (tronc), 2 (cuisse) ou 3(tibia)
%pmasse donne le rapport entre la masse segmentaire et la masse totale
pmasse = [0.497,0.1,0.0465];
%coef donne la position du CM segmentaire par rapport au point proximal sur la longueur du segment
coef = [0.495,0.433,0.433];
%rayon des segments corporels sur la hauteur
rg = [0.144,0.121,0.114];
%longueur des segments sur la hauteur totale
long = [0.288,0.245,0.246] ;

%calcul de la masse du segment
m = mglobal * pmasse(membres);
% calcul de la longueur du segment
l = hauteur * long(membres);
%rayon du segment corporel
r = rg(membres)*hauteur;

% matrice d'inertie passant par le milieu du cylindre m:
ineri1 = 0.25.*m.*r.^2+(1/12).*m.*l.^2 ;
ineri3 = 0.5.*m.*r.^2 ;
%utilisation du theoreme de huyguens pour le calcule de l'inertie (theoreme des axes paralleles
iner1 = ineri1 +m*(l*(0.5-coef(membres)))^2;
%valeur de laz matrice d'inertie
Ii=[iner1 0 0;0 iner1 0;0 0 ineri3];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FONCTION MASSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = masse(mglobal,membres,pmasse)

% membres = 1 (tronc), 2 (cuisse) ou 3(tibia)
%pmasse donne le rapport entre la masse segmentaire et la masse totale

%calcul de la masse du segment
m = mglobal * pmasse(membres);

