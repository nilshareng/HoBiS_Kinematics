%%% Le gros morceau...

% Déclaration de variables / flags
VarDeclaration;

dp = [0.001, 0.001];

% Ajout nécessaire, c'était pour ne pas modifier la fonction de calcul de
% la Jacobienne...
dp=[dp,Period];

% Affichages si désirés...
printflag=0;

% Ajout de l'option gradient 0, non updaté pour les c3d actuels car pas de convergence de base.
derflag = 0;
MaxLoop = 100;
SaveV = [];
SaveM = [];
% close all;
while c<MaxLoop % Nombre de cycle arbitraire, 20-25 suffisant pour discerner les convergences sans crash
    % Initialisation des Polynômes des TA.
    NPolA = [];
     
    % Update des PCA en ajoutant Modifs - Initialisé à 0
    for i = 1:size(PCA,1)
        NPCA(i,2:3) = NPCA(i,2:3) + [Modifs(2*i-1) , Modifs(2*i)];
        if c>1 % Ajout de l'erreur commise en cas de PCA éloignés de force 
            NPCA(i,2:3) = NPCA(i,2:3) + [InvE(i,c) ,0];
        end
        % 1-périodicité des PCA sur t
        NPCA(i,2) = NPCA(i,2) - fix(NPCA(i,2));
    end
    
    % Passage sur un intervalle positif en vue de retourner sur [0 1]
    NPCA(NPCA(:,2)<0,2) = NPCA(NPCA(:,2)<0,2)+1;
    
    % Compteur parce que la forme des PCA c'est dla merde pour certains calculs
    count=0;
       
    % Si on a touché aux PCA, on risque 2 PCA trop proches
    if c>1
        for i = 1:s
            % Pour chaque DDL qui n'est pas symétrisé
            tmP = NPCA(NPCA(:,1)==i,:);
            % On prend tous les PC d'un DDL, et on complère l'intervalle avec le premier +1 en t 
            tmP = [tmP;tmP(1,1),tmP(1,2)+1,tmP(1,3)];
            for j=1:size(tmP,1)-1
                % Pour les PC de ce DDL
                count = count +1;
                % Si le PC courant et son voisin suivant  sont trop proches
                if(abs(tmP(j,2) - tmP(j+1,2))<Threshold)
                    % Stockage de l'occurence, du cycle et du PC
                    Iflag = [Iflag, [c;i;j]];
                    if(tmP(j,2)-tmP(j+1,2)>0)
                        % S'ils ont été inversés lors des Modifs, remise dans l'ordre 
                        a=[tmP(j+1,:);tmP(j,:)];
                        tmP(j,:)= a(1,:);
                        tmP(j+1,:)=a(2,:);
                    end
                    
                    
                    if (j==size(tmP,1)-1) % Exception dans l'ordre s'il s'agit du dernier, vu que c'est le miroir du premier 
                        % Décalage des PCA trop proches
                        tmP(1,2)=tmP(1,2)+Threshold - (tmP(end,2) - tmP(j,2));
                        % Stockage et cumul des erreurs commises
                        InvE(count-j+1,c+1)=InvE(count-j+1,c)+abs(tmP(j,2) - tmP(j+1,2));
                    else % idem mais dans le cas général
                        tmP(j+1,2)=tmP(j+1,2)+Threshold - (tmP(j+1,2) - tmP(j,2));
                        InvE(count+1,c+1)=InvE(count+1,c)+abs(tmP(j,2) - tmP(j+1,2));
                    end
                end
            end
            % On récupère les PCA ordonnés, et séparés par au moins Threshold sur t 
            NPCA(NPCA(:,1)==i,:) = tmP(1:end-1,:);
        end
    end
    
    % PCA sur DDL 1 à 7 Séparés et ordonnées
    
    
    % Symétrisation Hard des PCA
    for i=8:11
        % On symétrise en déphasant
        NPCA(NPCA(:,1)==i,2:3)=[NPCA(NPCA(:,1)==i-4,2)+mid/Period, NPCA(NPCA(:,1)==i-4,3)];
        if i==9 || i==8
            % Et en inversant pour HY et HX
            NPCA(NPCA(:,1)==i,3)= -NPCA(NPCA(:,1)==i,3);
        end
        % Rangement des PCA dans l'ordre croissant de t
        NPCA(NPCA(:,1)==i,2) = NPCA(NPCA(:,1)==i,2)- fix(NPCA(NPCA(:,1)==i,2));
        tmp = NPCA(NPCA(:,1)==i,2:3);
        [~,I]=sort(tmp(:,2));
        NPCA(NPCA(:,1)==i,2:3) = tmp(I,:);
    end
 
    % Détection et rangement en cas d'inversion n°2
    for i=1:11
        tmp = NPCA(NPCA(:,1)==i,:);
        [tmp2,I] = sort(tmp(:,2));
        if (tmp2~=tmp(:,2))
            Cflag = [Cflag,[c;i]];
        end
        NPCA(NPCA(:,1)==i,:)=tmp(I,:);
    end
    
    % Stockage du NPCA calculé à ce cycle
    Storing = [Storing ; c*ones(size(NPCA,1),1) , NPCA];
    
    
% %  Gestion des détails - Deprecated  
%     if c>=0
%         tem = Details2(TempPCA, NPCA, E);
%         if isempty(tem) && ~isempty(Emem)
%             E = Emem(2:end, end-size(EpsilonOrig,2):end);
%         elseif isempty(tem)
%             E = EpsilonOrig;
%         else
%             E = tem;
%         end
%         Emem = [Emem , [ones(1,size(E,2));E]];
%     end

    % To prevent crashes... Toujours le pb de t =0 pour un PCA...
    for i =1:size(NPCA,1)
        if NPCA(i,2)==0
            NPCA(i,2)=eps;
        end
    end
    
    
    
    
    % Calcul des Polynomes NPolA entre les NPCA
    for i = 1:11
        temp = PC_to_spline(NPCA(NPCA(:,1)==i,2:3),1);
        NPolA = [NPolA ; i*ones(size(temp,1),1), temp(:,2:end)];
    end
    
    % Calcul de la fonction de coût pour les NPCA
    NRef = [];
    if derflag 
        DRef=[];
        [a,b] = ErrorFun4(NPolA,X,Sequence,Rmarkers,RReperes); % ErrorFun pour le gradient
        
        % Mesure et mémorisation de la norme de la somme des erreurs
        for i =1:size(X,1)
            NRef = [NRef; [a(i,:), b(i)]];
        end
    else
        [a,NewMarkers,NewReperes] = ErrorFun3(NPolA,X,Sequence,Rmarkers,RReperes); % ErrorFun2 sinon
        
        % Mesure et mémorisation de la norme de la somme des erreurs
        for i =1:size(X,1)
            NRef = [NRef; X(i,3:end)-a(i,:)];
        end
        
    end
    
    % Echantillonage pour calcul coûts secondaires
    [P, TA] = Sampling_txt(NPolA,Period,Sequence,Markers,Reperes);
    
    % Sécurité
    if(isempty(P)||isempty(TA))
        return;
    end
    
    % Affichages si printflag
    if (printflag)
        figure(1); % TA - Original, puis successifs
        hold on;
        for i = 1:11
            subplot(3,4,i);
            hold on;
            xlabel('cycle de marche, %');
            ylabel('trajectoire articulaire, rad');
            plot(0:1/Period:1,TA(:,i),'b');
            plot(0:1/Period:1,TrajAng(:,i),'b--');
        end
        sgtitle('Traj. angulaires sur cycle de marche. G à D : bX, bY,bZ, hgX, hgY, hgZ, ggX, hdX, hdY, hdZ, gdX');
        hold off;
        
        figure(2); % Poulaines - Originale en rouge, puis successifs, et X
        for i = 1:6
            subplot(2,3,i);
            hold on;
            xlabel('cycle de marche, %');
            ylabel('position cheville, m');
%             plot(0:1/(Period-1):1,Spline(:,i),'g--');
            plot(0:1/(Period-1):1,GT(:,i),'r');
%             plot(0:1/(Period-1):1,PN(:,i),'b');
            
            if i<4
                plot(X(1:size(X,1)/2,2),X(1:size(X,1)/2,i+2),'rx','MarkerSize',10);
                
            else
                plot(X(size(X,1)/2+1:end,2),X(size(X,1)/2+1:end,i-1),'rx','MarkerSize',10);
            end
            if c==0 && size(dp,2)==3
%                 plot(0:1/Period:1,PN(:,i));
                plot(0:1/(Period-1):1,P(1:end-1,i));
            else
                plot(0:1/(Period-1):1,P(1:end-1,i));
            end
            
        end
        sgtitle('Poulaines respectivement Gauche/Droite, X, Y, Z');
%         legend(leg)
    end
    
    % Coût énergétique
%     Cost = 0;
    Cost = ECShort(P,TA,M,Markers,Inertie);
    
    % Jerk
    JerkRef = 0;
%     Jerk(NPolA);
    
    % Sauvegarde des valeurs de contrôle
    Conv = [Conv , [norm(NRef); Cost ; JerkRef]];
    % Affichage
    if printflag
        norm(NRef)
        Cost
        JerkRef
    end
        
    
    % Calcul de la Jacobienne
    Manipflag =0;
    if Manipflag 
        Manip;
    end
%     [Jk, V, DJerk] = calc_jacobien_PC_3D(NPCA, NPolA,derflag, 0, X, dp,M, Cost, JerkRef, Fem1g, Fem1d, Fem6g, Fem6d, Tal1g, Tal1d, R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local);
    [Jk, V, DJerk] = calc_jacobien_PC_4D(NPCA, NPolA,X,dp,M,Cost,JerkRef,Rmarkers,RReperes,Sequence,Inertie);
    SaveV = [SaveV , V];
    SaveM = [SaveM ; NPolA];
    % Si problème de taille de tableau, complétion avec des 1
    if size(V,1)~=2*size(NPCA,1)
        V = [V ; ones(size(Modifs,1)-size(V,1),1)];
        DJerk = [DJerk ; ones(size(Modifs,1)-size(DJerk,1),1)];
    end
    
    % Si Erreur -> Exit
    if(isempty(Jk)||isempty(V)||isempty(DJerk))
%         break
    error('Jacobian is fucked');
    end
    
    % Pseudo Inverse
    Jkp = pinv(Jk);
    
    % Projection
    Proj = eye(2*size(PCA,1)) - Jkp*Jk;
    
    
    % Mémorisation
    N = [N , norm(NRef)];
    mV = [mV , Cost];
    mJ=[mJ, JerkRef];
    pModifs = [pModifs , Modifs];
    
    
    % Adaptation de la forme du vecteur de distance à X
    if derflag % Pas utilisé
        t=zeros(4*size(NRef,1),1);
        for a = 1:size(NRef,1)
            t(4*a-3)= X(a,3) - NRef(a,1);
            t(4*a-2)= X(a,4) - NRef(a,2);
            t(4*a-1)= X(a,5) - NRef(a,3);
            t(4*a)= - NRef(a,4);
        end
        delta = t;
        
    else % Cas courant
        t=zeros(3*size(NRef,1),1);
        for a = 1:size(NRef,1)
            t(3*a-2)=NRef(a,1);
            t(3*a-1)=NRef(a,2);
            t(3*a)=NRef(a,3);
        end
        delta = t;
    end
    
    delta = (delta/norm(delta))*min(norm(delta),0.01);
    
    
    % Coefficient pour pondérer les tâches secondaires -> Meme norme que
    % Proj
    weightCE = ceil(log(mean(mean(abs(Proj)))/(mean(abs(V)))));
    weightJerk = ceil(log(mean(mean(abs(Proj)))/(mean(abs(DJerk)))));
   
    % Si cassé
    if (isnan(weightCE)||isnan(weightJerk))  % Changements trop mineurs et Tâches secondaires Infinies par endroits
        Modifs = Jkp*delta;
    else
        
        % Pondération par rapport au log de la norme du vecteur de
        % modifications
        weightKer = ceil(log10(mean(mean(abs(Jkp*delta)))));
        if derflag
            Modifs = Jkp*delta;  %Proj*(VSter*(10^(weightKer+weightCE))*V + JSter*(10^(weightKer+weightJerk))*DJerk);
        else
            % En général, VSter et JSter =0, donc pas d'influence
%             Modifs = Jkp*delta + Proj*(V/norm(V))*VSter;
            VSter = 0.5*norm(Jkp*delta);
            Modifs = Jkp*delta + Proj*(V/norm(V))*VSter;
        end
    end
    % Sauvegarde des poids -> inutile ici car je ne les fais plus varier dynamiquement actuellement 
    mSter=[mSter ; VSter,JSter];
       
    % Affichage de la convergence
    if printflag
        figure(3);
        hold on;
        subplot(3,1,1);
        plot(1:size(Conv,2),Conv(1,1:end),'rx');
        subplot(3,1,2);
        plot(1:size(Conv,2),Conv(2,1:end),'rx');
        subplot(3,1,3);
        plot(1:size(Conv,2),Conv(3,1:end),'rx');
    end
    % Incrémentation du compteur de cycles
    c = c+1;
    fConv = figure(12);
    hold on;
    title('Convergence de la boucle d''optimisation - somme des distances aux empreintes - en m');
    subplot(2,1,1);
    hold on;
    plot(1:size(Conv,2),Conv(1,:));
    subplot(2,1,2);
    plot(1:size(Conv,2),Conv(2,:));
end

[Gait, GaitMarkers, GaitReperes] = Angles2Gait(TA,Sequence,Markers,Reperes,model.gait*1000,P*1000,SplinedComputedPoulaine*1000,X(:,3:5)*1000);%[model.gait(4:6,:)*1000, model.gait(1:3,:)*1000]);
fFGait = DisplayGait(GaitMarkers);

Results = struct;
Results.Inputs = answer;
Results.InitialPoulaine = model.gait;
Results.InitialReference = model.reference;
Results.InitialDescription = model.description;
Results.InitialSplinedPoulaine = SplinedComputedPoulaine;
Results.InitialSplinedAngles = SplinedAngles;
Results.X = X;
Results.InitialPolynom = PolA;
Results.FinalPolynom = NPolA;
Results.IncrementalPCModification = pModifs;
Results.Convergence = Conv;
Results.NCycles = c;
Results.Figure.Conv = fConv;
Results.Figure.FinalGait = fFGait;
Results.MemoryEC = SaveV;
Results.MemoryPoulaine = SaveM;


DisplayGait(GaitMarkers,20,'1');

% NPolA = SaveM(287:287+25,:);
% [P, TA] = Sampling_txt(NPolA,Period,Sequence,Markers,Reperes);
% [Gait, GaitMarkers, GaitReperes] = Angles2Gait(TA,Sequence,Markers,Reperes,model.gait*1000,P*1000,SplinedComputedPoulaine*1000,X(:,3:5)*1000);%[model.gait(4:6,:)*1000, model.gait(1:3,:)*1000]);
% fFGait = DisplayGait(GaitMarkers);
% DisplayCurves(P);
% DisplayCurves(TA);

%%
%%% Selecting, Culling, Saving, ...
% PostProcessing;


