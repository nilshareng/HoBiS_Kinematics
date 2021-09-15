    %%% MàJ 06/09/2021
%%% Script principal pour le lancer de batchs de tests
%%% Entrées : Chemins précisés ci-dessous (p et PathPreSet) menant aux données c3d et aux trajectoires cibles
%%% Sorties : Données issues de la boucle d'optimisation, sauvegardées en .mat dans SavePath
clear all;
close all;
clc;

%addpath("C:\Users\nhareng\Desktop\CodeCommente\hobis");

flag = struct;

flag.c3d = 0;
flag.txt = 0;
flag.dyn = 0;
flag.prints =0;


I = input("Prompt path to data directory and press Enter, default (press Enter) is C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\BDD\ \n", 's');

if isempty(I)
    % Chemin d'accès par défaut : .c3d et à l'excel compilant les données de marche :  
    p = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\BDD\';
    
    % Récupération des données de l'excel: 
    % Format : Nom fichier / Frame début et fin cycle de marche / Frame initiale
    % BornesMarches = Chiffres (BorneSup,BorneInf,FrameIni) ; Names = Noms
    [BornesMarche, Names] = xlsread(strcat(p,'Classement_Pas.xlsx'),'A2:D79');
    flag.c3d = 1;
else
    % Chemin d'accès spécifié par l'utilisateur
    p = I;
    p='C:\Users\nhareng\Documents\Prog\Matlab\Test\warping 2006\';
    a = input('Prompt file aggregator Name (.xlsx) - or simply single file name (.txt)\n','s');
    switch a(end-2:end)
        case 'lsx'
            [Names, Footprints]= xlsread(strcat(p,'a'),'A2:B79');
            flag.c3d = 1;
        case 'txt'
            Data = load(strcat(p,a));
            Names ={a(1:end-4)};
            flag.txt = 1;
    end       
end


% Déf du dossier de récéption des données calculées pour ce batch
I = input("Prompt path to save directory, default (press Enter) is C:\Users\nhareng\Desktop\CodeCommente\hobis\Resultats\Batch\1\ \n", 's');

if isempty(I)
    SavePath = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Resultats\Batch\1\';
else
    SavePath = I;
end

I = input("Prompt path to directory, default (press Enter) is C:\Users\nhareng\Desktop\CodeCommente\hobis\Resultats\Batch\NewPresets\ \n", 's');

if isempty(I)
    % Déf du dossier contenant les données précalculées :
    % i.e. les Poulaines utilisées comme cible, elles mêmes issues des
    % trajectoires angulaires de chaque fichier .c3d de la BDD
    % TL&DR : Les trajectoires initiales déformées par l'algo et ciblées
    % dans certains cas
    PathPreSet = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Resultats\Batch\NewPresets\';
    flag.presets = 0;
elseif I=='0'
    flag.presets = 1;
else
    PathPreSet = I;
    flag.presets = 0;    
end

% Si pas de données précalculées, faire tourner une fois en mettant presets
% à 1. - Dégueu, mais utilisé pour l'instant 

% Première phase : Pour chaque fichier .c3d répertorié dans l'excel : -> Extraire les
% données des marqueurs pendant le cycle de marche, lancer l'algo, puis sauver le résultat.
for ii=1:length(Names)%1%:11:length(Names)
    
    % Création d'un dossier pour cette marche - là où vont être stockées toutes les données calculées
    dirname= strcat(SavePath,Names{ii});
    mkdir(dirname);
    
    
    if flag.c3d
        
        % Bornes = Frames Début/Fin cycle de marche. Différence avec la frame initiale pour l'initialisation
        Bornes=[BornesMarche(ii,1)-BornesMarche(ii,3) BornesMarche(ii,2)-BornesMarche(ii,3)];
        
        % Taille de l'échantillon de frames à sélectionner
        Period = Bornes(2)-Bornes(1)+1;
        
        % Structure contenant les positions de marqueurs
        markers=struct;
        
        % Acquisition des marqueurs à partir des .c3d en 2étapes
        % 1 - Acquisition
        TN{ii}=strcat(Names{ii},'.c3d');
        NewFile = strcat(p,TN{ii});
        a = btkReadAcquisition(NewFile);
        premarkers = btkGetMarkers(a);
        Tmp = struct2cell(premarkers);
        
        % 2 - Etape intermédiaire : Les noms, ordres, et nombres des marqueurs varient à chaque .c3d
        % Il faut donc repérer et sélectionner ceux qui nous intéresse.
        % Sélection G/D - FWT BWT KNE KNI ANE ANI
        
        % Récupération des noms des marqueurs
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
                    case 'REMAL' % REMAL et RIMAL inversés pour normal_walking_4 
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
        % Marqueurs acquis, Lancement de la première phase de l'algorithme :
        % Créations des repères articulaires, Filtrage, Cinématique Inverse et
        % Approximation par Splines
        
        C3DUpdate2;
        
    elseif flag.txt
        Period = 100;
        [markers, Param] = HobisDataParser(Data);
        C3DUpdate2;
    end
    
    if flag.presets
        close all;
        save(strcat(PathPreSet,Names{ii},'.mat'),'PN','Pol','Param', 'R_monde_local','R_Pelvis_monde_local', 'R_LFem_ref_local', 'R_LTib_ref_local', 'R_RFem_ref_local', 'R_RTib_ref_local');
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
        for jj =1:length(Names)%1:length(Names)-1
            jj
            % Deuxième phase, à la marche courante, on va successivement appliquer
            % les autres marches en tant que cibles.
            if size(Names{jj}(1:end-3),2)==size(Names{ii}(1:end-3),2) || flag.dyn
                if flag.dyn %Names{jj}(1:end-3)== Names{ii}(1:end-3) || flag.dyn
                    PreLoop;
                    
                    % Troisième phase, Boucle d'optimisation.
                    Loop_Batch;
                    close all;
                    
                    % Sauvegarde des données dans des .mat importables
                    save(strcat(dirname,'\',Names{jj},'P','.mat'),'Param','Saved','SNPCA','GT','Conv','PFin','TAFin','X','mem','Iflag','Storing')
                else
                    PreLoop;
                    Loop_Batch;
                    close all; % Just in case
                    
                    save(strcat(dirname,'\',Names{jj},'P','.mat'),'Param','Saved','SNPCA','GT','Conv','PFin','TAFin','X','mem','Iflag','Storing')
                end
            end
        end
    end
%     % Si besoin ici, sauvegarde des PreSets, ie données utilisées comme cibles
%     % dans la phase 2 et 3

end
    
    
    