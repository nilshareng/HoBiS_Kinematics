%%% M�J 06/09/2021
%%% Script principal pour le lancer de batchs de tests
%%% Entr�es : Chemins pr�cis�s ci-dessous (p et PathPreSet) menant aux donn�es c3d et aux trajectoires cibles
%%% Sorties : Donn�es issues de la boucle d'optimisation, sauvegard�es en .mat dans SavePath
clear all;
close all;
clc;

%addpath("C:\Users\nhareng\Desktop\CodeCommente\hobis");

flag = struct;

flag.c3d = 0;
flag.txt = 0;
flag.dyn = 0;
flag.prints =0;

%% Input part

prompt = ...
    {'Enter Ressources directory (absolute) path:',...
    'Enter Description Posture file (absolute) path:',...
    'Enter Reference Posture file (absolute) path:',...
    'Enter PreSets files directory (absolute) path:',...
    'Initial Poulaine (absolute) path:',...
    'Model Mass:', ...
    'Model Articular Range of Motion Min Boundaries in degrees (Pelv*3, LHip*3, LKnee*1, RHip*3, RKnee*1):',...
    'Model Articular Range of Motion Max Boundaries in degrees (Pelv*3, LHip*3, LKnee*1, RHip*3, RKnee*1):',...
    'Footprints (target) left footstrike (gait% X Y Z):',...
    'Footprints (target) left toe off (gait% X Y Z):',...
    'Footprints (target) left min Z (gait% X Y Z):'...
    'Footprints (target) right footstrike (gait% X Y Z):',...
    'Footprints (target) right toe off (gait% X Y Z):',...
    'Footprints (target) right min Z (gait% X Y Z):'...
    'Enter Save folder Path:'};
dlgtitle = 'Inputs';
dims = [1 110];
definput = {'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\description.txt',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\reference.txt',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\NewestPresets\',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\poulaineAppui.txt',...
    '70',...
    '[45 45 45 20 45 30 110 20 45 30 110]',...
    '[-45 -45 -45 -90 -30 -60 15 -90 -30 -60 15]',...
    '[0.2000    0.0650   -0.0188   -0.7755]',...
    '[0.7333    0.1135    0.0585   -0.6367]',...
    '[0.2833    0.0656    0.0595   -0.7840]'...
    '[0.7333   -0.0483   -0.0049   -0.7827]',...
    '[0.2667   -0.0990    0.0351   -0.6603]',...
    '[0.8167   -0.0545    0.0777   -0.7835]'...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Resultats\Txt'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

model = struct;

p = answer{1};
DataDes = load(answer{2});
model.description = DataDes;

DataRef = load(answer{3});
model.reference = DataRef;

if isempty(answer{4})
    flag.presets = 1;
    mkdir p 'Presets'
    PathPreSet = strcat(answer{1},'Presets\');
else
    PathPreSet = answer{4};
    flag.presets = 0;
end

model.gait = load(answer{5});

M = str2double(answer{6});

model.jointRangesMin = str2num(answer{7});
model.jointRangesMax = str2num(answer{8});

X = [str2num(answer{9}) ; str2num(answer{10}) ; str2num(answer{11}); str2num(answer{12}) ; str2num(answer{13}) ; str2num(answer{14})];
X = [[2;2;3;5;5;6] , X];
flag.steps = 1;
flag.txt = 1;
SavePath = answer{15};
Names = {DataDes(1:end-4)};


%%% Variables to set : p , PathPreSet, Desc


%%
% 
% 
% I = input("Prompt path to data directory and press Enter, default (press Enter) is C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\BDD\ \n", 's');
% 
% if isempty(I) %|| 1
%     % Chemin d'acc�s par d�faut : .c3d et � l'excel compilant les donn�es de marche :
%     p = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\BDD\';
%     
%     % R�cup�ration des donn�es de l'excel:
%     % Format : Nom fichier / Frame d�but et fin cycle de marche / Frame initiale
%     % BornesMarches = Chiffres (BorneSup,BorneInf,FrameIni) ; Names = Noms
%     [BornesMarche, Names] = xlsread(strcat(p,'Classement_Pas.xlsx'),'A2:D79');
%     flag.c3d = 1;
% else
%     % Chemin d'acc�s sp�cifi� par l'utilisateur
%     p = I;
%     p='C:\Users\nhareng\Documents\Prog\Matlab\Test\warping 2006\';
%     p = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\';
%     a = input('Prompt Des file aggregator Name (.xlsx) - or simply single file name des (.txt)\n','s');
%     b = input('Prompt Ref file aggregator Name (.xlsx) - or simply single file name ref (.txt)\n','s');
%     switch a(end-2:end)
%         case 'lsx'
%             [Names, Footprints]= xlsread(strcat(p,a),'A2:B79');
%             flag.c3d = 1;
%         case 'txt'
%             DataDes = load(strcat(p,a));
%             DataRef = load(strcat(p,b));
%             Names ={a(1:end-4)};
%             flag.txt = 1;
%     end
% end
% %  p = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\BDD\';
% 
% % D�f du dossier de r�c�ption des donn�es calcul�es pour ce batch
% I = input("Prompt path to save directory, default (press Enter) is C:\Users\nhareng\Desktop\CodeCommente\hobis\Resultats\Batch\1\ \n", 's');
% 
% if isempty(I)
%     SavePath = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Resultats\Batch\1\';
% else
%     SavePath = I;
% end
% 
% %%% Need to setup 'p' 'SavePath' 'flag.txt' 'PathPreSet' 'flag.presets'
% 
% 
% I = input("Prompt path to presets directory, default (press Enter) is C:\Users\nhareng\Desktop\CodeCommente\hobis\Resultats\Batch\NewPresets\ \n", 's');
% 
% if isempty(I)
%     % D�f du dossier contenant les donn�es pr�calcul�es :
%     % i.e. les Poulaines utilis�es comme cible, elles m�mes issues des
%     % trajectoires angulaires de chaque fichier .c3d de la BDD
%     % TL&DR : Les trajectoires initiales d�form�es par l'algo et cibl�es
%     % dans certains cas
%     PathPreSet = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\NewestPresets\';
%     flag.presets = 0;
% elseif I=='0'
%     PathPreSet = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\NewestPresets\';
%     flag.presets = 1;
% else
%     PathPreSet = I;
%     flag.presets = 0;
% end

%%

% Si pas de donn�es pr�calcul�es, faire tourner une fois en mettant presets
% � 1. - D�gueu, mais utilis� pour l'instant

% Premi�re phase : Pour chaque fichier .c3d r�pertori� dans l'excel : -> Extraire les
% donn�es des marqueurs pendant le cycle de marche, lancer l'algo, puis sauver le r�sultat.
for ii=1:length(Names)%1%:11:length(Names)
    
    % Cr�ation d'un dossier pour cette marche - l� o� vont �tre stock�es toutes les donn�es calcul�es
    dirname= strcat(SavePath,Names{ii});
%     mkdir(dirname);
    
    
    if flag.c3d
        
        % Bornes = Frames D�but/Fin cycle de marche. Diff�rence avec la frame initiale pour l'initialisation
        Bornes=[BornesMarche(ii,1)-BornesMarche(ii,3) BornesMarche(ii,2)-BornesMarche(ii,3)];
        
        % Taille de l'�chantillon de frames � s�lectionner
        Period = Bornes(2)-Bornes(1)+1;
        
        % Structure contenant les positions de marqueurs
        markers=struct;
        
        % Acquisition des marqueurs � partir des .c3d en 2�tapes
        % 1 - Acquisition
        TN{ii}=strcat(Names{ii},'.c3d');
        NewFile = strcat(p,TN{ii});
        a = btkReadAcquisition(NewFile);
        premarkers = btkGetMarkers(a);
        Tmp = struct2cell(premarkers);
        
        % 2 - Etape interm�diaire : Les noms, ordres, et nombres des marqueurs varient � chaque .c3d
        % Il faut donc rep�rer et s�lectionner ceux qui nous int�resse.
        % S�lection G/D - FWT BWT KNE KNI ANE ANI
        
        % R�cup�ration des noms des marqueurs
        b = fieldnames(premarkers);
        if flag.dyn
            for j = 1:size(b,1)
                tempName = b{j};
                if size(tempName,2) > 4
                    tempName = tempName(end-3:end);
                end
                %             tempName = tempName(end-3:end);
                switch tempName
                    case 'RASIS'
                        markers.RFWT = Tmp{j};
                    case 'LASIS'
                        markers.LFWT = Tmp{j};
                    case 'RPSIS'
                        markers.RBWT = Tmp{j};
                    case 'LPSIS'
                        markers.LBWT = Tmp{j};
                    case 'LECE'
                        markers.LKNE = Tmp{j};
                    case 'RECE'
                        markers.RKNE = Tmp{j};
                    case 'LICE'
                        markers.LKNI = Tmp{j};
                    case 'RICE'
                        markers.RKNI = Tmp{j};
                    case 'LEMAL'
                        markers.LANE = Tmp{j};
                    case 'RIMAL'
                        markers.RANE = Tmp{j};
                    case 'LIMAL'
                        markers.LANI = Tmp{j};
                    case 'REMAL' % REMAL et RIMAL invers�s pour normal_walking_4
                        markers.RANI = Tmp{j};
                        
                        
                end
            end
            
        elseif flag.c3d
            for j = 3:size(b,1)
                tempName = b{j};
                if size(tempName,2) > 4
                    tempName = tempName(end-3:end);
                end
                switch tempName
                    case 'RFWT'
                        markers.RFWT = Tmp{j};
                    case 'LFWT'
                        markers.LFWT = Tmp{j};
                    case 'RBWT'
                        markers.RBWT = Tmp{j};
                    case 'LBWT'
                        markers.LBWT = Tmp{j};
                    case 'LKNE'
                        markers.LKNE = Tmp{j};
                    case 'RKNE'
                        markers.RKNE = Tmp{j};
                    case 'LKNI'
                        markers.LKNI = Tmp{j};
                    case 'RKNI'
                        markers.RKNI = Tmp{j};
                    case 'LANE'
                        markers.LANE = Tmp{j};
                    case 'RANE'
                        markers.RANE = Tmp{j};
                    case 'LANI'
                        markers.LANI = Tmp{j};
                    case 'RANI'
                        markers.RANI = Tmp{j};
                        
                        
                end
                
                
                
            end
            
        end
        % Marqueurs acquis, Lancement de la premi�re phase de l'algorithme :
        % Cr�ations des rep�res articulaires, Filtrage, Cin�matique Inverse et
        % Approximation par Splines
        
        C3DUpdate2;
        
    elseif flag.txt
        Period = 100;
%         [Dmarkers, DParam ]= HobisDataParser(DataDes);
        [Rmarkers, RParam ] = HobisDataParser(DataRef);
%         [DReperes DSeq NDmarkers NDParam] = ReperesFromMarkers(Dmarkers);
        [RReperes, RSeq, NRmarkers, NRParam] = ReperesFromMarkersCorrected(Rmarkers);
        A = zeros(1,11);
%         [tmp, tmpS, tmpmarkers] = fcineshort(A,NDParam, DReperes, Dmarkers);
%         [tmp, tmpS, tmpmarkers] = fcineshort(A,NDParam, DReperes, NDmarkers);
        
%         ConditionA0 = all(all(round(NRParam' - tmpS(2:end,:),6)==zeros(size(NRParam'))));
%         if ~ConditionA0
%         end
%         tmp*1000 - [NDmarkers.LTal1'; NDmarkers.RTal1']
%         [tmp, tmpS, tmpmarkers] = fcineshort(A,NRParam, RReperes, Rmarkers);
%         [tmp, tmpS, tmpmarkers] = fcineshort(A,NRParam, RReperes, NRmarkers);
        
        
        Sequence.Pelvis = 'xyz';
        Sequence.LHip = 'zyx';
        Sequence.LKnee = 'z';
        Sequence.RHip = 'zyx';
        Sequence.RKnee = 'z';
        
        [PosC,Markers,Reperes] = fcinematique([0 0 0 0 0 0 0 0 0 0 0], Sequence, Rmarkers, RReperes);
        
%         c.LFemur1(1:3,4) = RReperes.LFemur1(1:3,4);
%         c.LFemur2(1:3,4) = RReperes.LFemur2(1:3,4);
%         c.RFemur2(1:3,4) = RReperes.RFemur2(1:3,4);
%         c.RFemur1(1:3,4) = RReperes.RFemur1(1:3,4);
%         c.LTibia(1:3,4) = RReperes.LTibia(1:3,4);
%         c.RTibia(1:3,4) = RReperes.RTibia(1:3,4);
        
        DisplayMarkers(Markers,1,Reperes);
        
%         ConditionA0 = all(all(round(NRParam' - tmpS(2:end,:),6)==zeros(size(NRParam'))));
%         if ~ConditionA0
%             error('Shit not lining up in kinematic function, angles0 =/= Initial posture')
%         end
%         tmp*1000 - [NRmarkers.LTal1'; NRmarkers.RTal1']
%         for i =1:11
%             A = zeros(1,11);
%             A(i)=pi/3
%             [tmp, tmpmarkers] = fcinematique(A, Sequence, Rmarkers, RReperes);
%             %                 [tmp, tmpS, tmpmarkers] = fcineshort(A,NDParam, DReperes, NDmarkers);
%             tmp*1000 - [NDmarkers.LTal1'; NDmarkers.RTal1']
%             %                 figure(i);
%             DisplayMarkers(tmpmarkers);
%             A = zeros(1,11);
%             A(i)=2*pi/3
%             [tmp, tmpmarkers] = fcinematique(A, Sequence, Rmarkers, RReperes);
%             tmp*1000 - [NDmarkers.LTal1'; NDmarkers.RTal1']
%             %                 figure(i);
%             DisplayMarkers(tmpmarkers);
%         end
        % TODO - Boucler sur Alterations, matcher la pos ini chev � partir
        % de poulaine charg�e et Posture Ref - ne partir sur pos descr que si bexoin compa
        % extra espece
%         Alterations;
%         close all;
%         PreLoop;
%         Loop_Batch_Txt;
%         close all;
%         save(strcat(dirname,'\',Names{jj},'P','.mat'),'Param','Saved','SNPCA','GT','Conv','PFin','TAFin','X','mem','Iflag','Storing')

    end
    if flag.presets
        close all;
        save(strcat(PathPreSet,Names{ii},'.mat'),'PN','Pol','Param', 'R_monde_local','R_Pelvis_monde_local', 'R_LFem_ref_local', 'R_LTib_ref_local', 'R_RFem_ref_local', 'R_RTib_ref_local');
        jj=ii;
        PreLoop;
        save(strcat(PathPreSet,Names{ii},'.mat'),'PN','Pol','Param','X', 'R_monde_local','R_Pelvis_monde_local', 'R_LFem_ref_local', 'R_LTib_ref_local', 'R_RFem_ref_local', 'R_RTib_ref_local');
        
        %     if flag.dyn
        %         Test = spline_to_curve(Pol, 1, 1/202 );
        %         Test = Test(1:end-3,:);
        %         Test =[Test(126:end,:);Test(1:125,:)];
        %         Test = [Test ; Test ; Test];
        %         fileID = fopen('C:\Users\nhareng\Desktop\CodeCommente\DataLouise\DataLouise.txt','w');
        %         fprintf(fileID,'%8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n', 'PelX', 'PelY', 'PelZ', 'RHipX', 'RHipY', 'RHipZ', 'LHipX', 'LHipY', 'LHipZ', 'RKneeX', 'LKneeX');
        %         fprintf(fileID,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n', Test);
        %         fclose(fileID);
        %     end
    else
        for jj =1%1:length(Names)%1:length(Names)-1
            %             jj
            % Deuxi�me phase, � la marche courante, on va successivement appliquer
            % les autres marches en tant que cibles.
            if size(Names{jj}(1:end-3),2)==size(Names{ii}(1:end-3),2) || flag.dyn
                if flag.dyn %Names{jj}(1:end-3)== Names{ii}(1:end-3) || flag.dyn
                    PreLoop;
                    
                    % Troisi�me phase, Boucle d'optimisation.
                    Loop_Batch;
                    close all;
                    
                    % Sauvegarde des donn�es dans des .mat importables
                    save(strcat(dirname,'\',Names{jj},'P','.mat'),'Param','Saved','SNPCA','GT','Conv','PFin','TAFin','X','mem','Iflag','Storing')
                else
                    if flag.txt
                        Alterations;
                        PreLoop;
                        Loop_Batch_Txt;
                        close all;
                    end
%                     PreLoop;
%                     Loop_Batch;
                    close all; % Just in case
                    
                    save(strcat(dirname,'\',Names{jj},'P','.mat'),'Param','Saved','SNPCA','GT','Conv','PFin','TAFin','X','mem','Iflag','Storing')
                end
            end
        end
    end
    %     % Si besoin ici, sauvegarde des PreSets, ie donn�es utilis�es comme cibles
    %     % dans la phase 2 et 3
    
end




