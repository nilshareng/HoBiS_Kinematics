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
SC = 0;
PasDelta = 0.1; % Norm of the vector from curr pos to target pos
PasModifs = 1; % 
close all;
fConv = figure(12);
hold on;

Mem = struct;

Log = struct;
Log.ModifsNorm = zeros(MaxLoop+1,1);
% TimeIniOpti = tic();

% Log.TimeIniCycle = zeros(MaxLoop+1,1);
% Log.TimeCycle = zeros(MaxLoop+1,1);

while c<MaxLoop % Nombre de cycle arbitraire, 20-25 suffisant pour discerner les convergences sans crash
    
%     Log.TimeIniCycle(c+1) = tic();
    
    c

    [NPCA,NPolA,Iflag] = ModifPCA(Modifs, NPCA, Period, Sequence, Markers, Reperes, InvE, mid,...
        Threshold, PCA, c, Cflag, rate, Iflag);

    [tmpP2, tmpTA2] = Sampling_txt(NPolA,Period,Sequence,Markers,Reperes);
    
    Mem(c+1).PC = NPCA;
    Mem(c+1).Pol = Pol;
    
%     DisplayCurves(tmpP2,10);
%     DisplayCurves(tmpTA2,11);
%     DisplayCurves(InitialGait,10);
%     DisplayX(X,size(tmpP2,1)-2,10)
%     
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
        % Cibles intermédiaires :
        [a,NewMarkers,NewReperes] = ErrorFun3(NPolA,X,Sequence,Rmarkers,RReperes); % ErrorFun2 sinon
        % Mesure et mémorisation de la norme de la somme des erreurs
        
        for i =1:size(X,1)
            NRef = [NRef; X(i,3:end)-a(i,:)];
        end
        
    end
    
    % Echantillonage pour calcul coûts secondaires
    [P, TA] = Sampling_txt(NPolA,Period,Sequence,Markers,Reperes);
    
    Mem(c+1).P = P;
    Mem(c+1).TA = TA;
    
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
            plot(0:1/(Period-1):1,TA(1:end-1,i),'b');
            plot(0:1/(Period-1):1,TrajAng(:,i),'r--');
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
    
    Mem(c+1).EC = Cost;
    
    % Jerk
    JerkRef = 0;
%     Jerk(NPolA);

    % ArticularCost 
    
    ACRef = ArticularCostPC(NPolA, Period, Sequence, Markers, Reperes, model.jointRangesMin, model.jointRangesMax);
    
    Mem(c+1).AC = ACRef;
    
    % Sauvegarde des valeurs de contrôle
    Conv = [Conv , [norm(NRef); Cost ; JerkRef ; norm(ACRef)]];
    % Affichage
    if printflag
        norm(NRef)
        Cost
        JerkRef
    end
       
    Mem(c+1).Conv = Conv;
    
    % Calcul de la Jacobienne
    Manipflag =0;
    if Manipflag 
        Manip;
    end
%     [Jk, V, DJerk] = calc_jacobien_PC_3D(NPCA, NPolA,derflag, 0, X, dp,M, Cost, JerkRef, Fem1g, Fem1d, Fem6g, Fem6d, Tal1g, Tal1d, R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local);
    [Jk, V, DJerk, DAC] = calc_jacobien_PC_4D(NPCA, NPolA,X,dp,M,Cost,JerkRef,ACRef,Rmarkers,RReperes, ...
        Sequence,Inertie, model.jointRangesMin, model.jointRangesMax);
    SaveV = [SaveV , V];
    SaveM = [SaveM ; NPolA];
    
    Mem(c+1).Jk = Jk;
    
    
    % Si problème de taille de tableau, complétion avec des 1
    if size(V,1)~=2*size(NPCA,1)
        V = [V ; ones(size(Modifs,1)-size(V,1),1)];
        DJerk = [DJerk ; ones(size(Modifs,1)-size(DJerk,1),1)];
    end
    
    % Si Erreur -> Exit
    if(isempty(Jk)||isempty(V)||isempty(DJerk))
%         break
    Log.FatalFunctionErrorCycle = c;
    error('Jacobian is fucked');
    end
    
    % Pseudo Inverse
    Jkp = pinv(Jk);
    
    
    % Protection contre singularités : 
    if norm(Jk) > 10^3 || norm(Jkp) > 10^3 
        [NPCA,NPolA,Jk, Jkp, V, DJerk, DAC, Log] = JacSingularitiesCheck( NPCA, NPolA, X,dp,M,Cost,JerkRef , ...
            ACRef, Period, Sequence, Markers, Reperes, InvE, mid, Threshold, PCA, c, Cflag, rate, ...
            model, tmpP2, tmpTA2,Log, SC, Jk, Jkp);
    end
    
    % Projection
    Proj = eye(2*size(PCA,1)) - Jkp*Jk;
    
    % Mémorisation
    N = [N , norm(NRef)];
    mV = [mV , Cost];
    mJ=[mJ, JerkRef];
    pModifs = [pModifs , Modifs];
    
    Mem(c+1).N = N;
    Mem(c+1).Modifs = Modifs;
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
    delta = (delta/norm(delta))*min(norm(delta),PasDelta);
    
%     for ij = 1:size(tmpTA2,1)
%         deltatheta2 = zeros(1,11);
%         tmp = tmpTA2(i,:) + ;
%         delta2 = fminsearch(@(deltatheta2) ArticularCost(tmp,deltatheta2,Proj,model.jointRangesMin, ...
%           model.jointRangesMax),deltatheta2);
%         tmp = tmp + (proj*delta2')';
%     end
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
            VSter = -0.2 * norm(Jkp*delta);
            PondAC = -0.5 * norm(Jkp*delta);
            
             %+ Proj*(V/norm(V))*VSter + Proj*(DAC/norm(DAC))*PondAC;
            TS1 = Proj * V / norm(Proj * V) * VSter;
            TS2 = Proj * DAC / norm(Proj * DAC) * PondAC;
            
            Modifs = Jkp*delta + TS1 ;%+ TS2;
            
%             %% 
%             
%             [TNPCA,TNPolA,Iflag] = ModifPCA(Modifs, NPCA, Period, Sequence, Markers, Reperes, InvE, ...
%                 mid, Threshold, PCA, c, Cflag, rate, Iflag);
%             
%             [Jk, V, DJerk, DAC] = calc_jacobien_PC_4D(TNPCA, TNPolA,X,dp,M,Cost,JerkRef,ACRef, ...
%                 Rmarkers,RReperes,Sequence,Inertie, model.jointRangesMin, model.jointRangesMax);
% 
%             Jkp = pinv(Jk);
%                         
%             % Protection contre singularités : 
%             if norm(Jk) > 10^3 || norm(Jkp) > 10^3
%                 [NPCA,NPolA,Jk, Jkp, V, DJerk, DAC, Log] = JacSingularitiesCheck(NPCA, NPolA, X,dp,M,Cost, ...
%                     JerkRef , ACRef, Period, Sequence, Markers, Reperes, InvE, mid, Threshold, PCA, c, Cflag, ...
%                     rate, model, tmpP2, tmpTA2,Log, SC, Jk, Jkp);
%             end
%             Fflag = 1;
%             if Fflag
%                 Proj = eye(2*size(PCA,1)) - Jkp*Jk;
%                 
%                 deltatheta2_Ini = 0.1 * ones(size(Modifs));
%                 options = optimset('Display','Iter','TolFun',1e-1,'MaxIter',5);
%                 deltatheta2 = fminsearch(@(deltatheta2) EnergyCostOptim(deltatheta2, Modifs, Proj, ...
%                     NPCA, Period, Sequence, Markers, Reperes, InvE, mid, Threshold, PCA, c, Cflag, ...
%                     rate, M, Inertie,Iflag), deltatheta2_Ini, options);
%                 
%                 Modifs = Modifs + (Proj*deltatheta2)*VSter;%Proj*(V/norm(V))*VSter + Proj*(DAC/norm(DAC))*PondAC;
%             end
%             %%
            
            
            

        end
    end
    Log.ModifsNorm(c+1) = norm(Modifs);
    Modifs = (Modifs/norm(Modifs))*min(norm(Modifs),PasModifs);
    
%     % Sauvegarde des poids -> inutile ici car je ne les fais plus varier dynamiquement actuellement 
%     mSter=[mSter ; VSter,JSter];
       
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
%     fConv = [];
%     fConv = figure(12);
%     hold on;
%     subplot(2,1,1);
%     hold on;
%     title('Convergence de la boucle d''optimisation - somme des distances aux empreintes - en m');
%     plot(1:size(Conv,2),Conv(1,:));
%     subplot(2,1,2);
%     hold on;
%     title('Convergence de la boucle d''optimisation - Travail des forces internes');
%     plot(1:size(Conv,2),Conv(2,:));
%     subplot(3,1,3);
%     hold on;
%     title('Convergence de la boucle d''optimisation - Cout Articulaire (C<10e20 pour rester dans les limites)');
%     plot(1:size(Conv,2),Conv(3,:));
    
    Log.OptiCycle.ClosureCycle = c;
%     Log.TimeCycle(c) = toc(Log.TimeIniCycle(c));
    
    if Conv(1,end) < 0.01
       c=MaxLoop+1; 
    else
        
    end
end

%% Selection de la meilleure solution

[~,d] = min(Conv(1,:));
Mem = Mem(d);

%%

% Log.TimeEndOpti = toc(Log.TimeIniOpti);
[Gait, GaitMarkers, GaitReperes] = Angles2Gait(tmpTA2,Sequence,Markers,Reperes,model.gait*1000, ...
    P*1000,SplinedComputedPoulaine*1000,X(:,3:5)*1000);%[model.gait(4:6,:)*1000, model.gait(1:3,:)*1000]);
% close(figure(1));
close all;
% fFGait = DisplayGait(GaitMarkers);

% DisplayCurves(P,2);       
% DisplayCurves([P(:,4:6), P(:,1:3)],2);
% DisplayCurves(tmpTA2,3);
% DisplayCurves([tmpTA2(:,1:3), tmpTA2(:,8:11) , tmpTA2(:,4:7)],3);

Results = struct;
Results.Inputs = answer;
Results.Mem = Mem;
Results.InitialPoulaine = InitialGait;
Results.InitialPoulaineScaled = model.gait;
Results.InitialPoulaineScaledIK = NewPoul;
Results.InitialSplinedPoulaine = SplinedComputedPoulaine;
Results.FinalPoulaine = P;
Results.TAPostIK = NewAngles;
Results.TAFinal = TA;
Results.InitialReference = model.reference;
Results.InitialDescription = model.description;
Results.InitialSplinedAngles = SplinedAngles;
Results.OriginalX = InputX;
Results.EmpreintesPoulaineInput = PX;
Results.PoulaineEmpreintesInput = OOPN;
Results.PoulaineEmpreintesInputScaled = OPN;
Results.EmpreintesScaled = X;
Results.InitialPolynom = PolA;
Results.FinalPolynom = NPolA;
Results.IncrementalPCModification = pModifs;
Results.Convergence = Conv;
Results.NCycles = c;
Results.Logs = Log;
if exist('PolX')
    Results.TargetPolynom = PolX;
%     Results.Target = ;
end

% Results.Figure.Conv = fConv;
% Results.Figure.FinalGait = fFGait;
Results.MemoryEC = mV;
Results.MemoryPoulaine = SaveM;
Results.GaitMarkers = GaitMarkers;

% DisplayGait(GaitMarkers,20,'4');

Spline = [];
I = 60;
for i = 1:11
    tmpPol = Pol(Pol(:,1)==i,:);
    TS2 = [];
    for tc= 0:1/I:1
        % Echantillonnage pour chaque Angle
        a= EvalSpline(tmpPol, tc);
        if(a==-10)
            P2=[];
            TA2=[];
           return; 
        end
        TS2 = [TS2;a];
    end
    Spline = [Spline, TS2];
end

