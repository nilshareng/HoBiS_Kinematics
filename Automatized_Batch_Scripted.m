
% Inertie.Rmasse = [0.497,0.1,0.0465];
% Inertie.CM = [0.495,0.433,0.433];
% Inertie.radius = [0.144,0.121,0.114];
% Inertie.length = [0.288,0.245,0.246] ;
% '[0.7000    -0.0650   -0.0188   -0.7755]',...
%     '[0.2333   -0.1135    0.0585   -0.6367]',...
%     '[0.7833   -0.0656    0.0595   -0.7840]'...
%     '[0.7333   -0.0483   -0.0049   -0.7827]',...
%     '[0.2667   -0.0990    0.0351   -0.6603]',...
%     '[0.8167   -0.0545    0.0777   -0.7835]'...

% 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\NewestPresets\antho056.mat',...
% 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\poulaine_base.txt',...
answer = definput;%inputdlg(prompt,dlgtitle,dims,definput);

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

InitialGaitPath = answer{5};
if size(InitialGaitPath,2)==7
    InitialGaitPath = InitialGaitPath(:,2:end);
end
% model.gait = InitialGaitPath;

M = str2double(answer{6});

model.jointRangesMin = str2num(answer{7})*pi/180;
model.jointRangesMax = str2num(answer{8})*pi/180;

% X = [str2num(answer{9}) ; str2num(answer{10}) ; str2num(answer{11})];%; str2num(answer{12}) ; str2num(answer{13}) ; str2num(answer{14})];
% X = [X ; X];
% X(4:6,:) = [X(4:6,1) + 0.5 , X(4:6,2)*-1 , X(4:6,3:4)];
% X(X(4:6,1)>1,1) = X(X(4:6,1)>1,1) -1;
% X = [[2;2;3;5;5;6] , X];
flag.steps = 1;
flag.txt = 1;
SavePath = answer{9};
Names = {DataDes(1:end-4)};

% prompt = {...
%     'Inertia - Pelvis segment, expressed at segment center of mass'...
%     'Inertia - Hip segment, expressed at segment center of mass'...
%     'Inertia - Leg segment, expressed at segment center of mass'...
%     'Inertia Center of mass position (% of segment - proximal to distal)'};
% 
% dlgtitle = 'Inertial inputs';
% dims = [1 110];
definput2 = {...
    '[60.8189 0 0 ; 0 60.8189 0 ; 0 0 52.1215]'...
    '[0.1807 0 0 ; 0 0.1807 0 ; 0 0  0.1481]'...
    '[0.0806 0 0 ; 0 0.0806 0 ; 0 0 0.0611]'...
    '[0.497,0.1,0.0465]'};

answerInertia = definput2;

Inertie.Pelvis = str2num(answerInertia{1});
Inertie.Hip = str2num(answerInertia{2});
Inertie.Leg = str2num(answerInertia{3});
Inertie.Coef = str2num(answerInertia{4});


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
        
        Sequence.Pelvis = 'xyz';
        Sequence.LHip = 'zyx';
        Sequence.LKnee = 'z';
        Sequence.RHip = 'zyx';
        Sequence.RKnee = 'z';
        
        [PosC,Markers,Reperes] = fcinematique([0 0 0 0 0 0 0 0 0 0 0], Sequence, Rmarkers, RReperes);
        
%         DisplayMarkers(Markers,1,Reperes);
        
    end
    if flag.presets
        close all;
        save(strcat(PathPreSet,Names{ii},'.mat'),'PN','Pol','Param', 'R_monde_local','R_Pelvis_monde_local', 'R_LFem_ref_local', 'R_LTib_ref_local', 'R_RFem_ref_local', 'R_RTib_ref_local');
        jj=ii;
        PreLoop;
        save(strcat(PathPreSet,Names{ii},'.mat'),'PN','Pol','Param','X', 'R_monde_local','R_Pelvis_monde_local', 'R_LFem_ref_local', 'R_LTib_ref_local', 'R_RFem_ref_local', 'R_RTib_ref_local');
    else
        Alterations;
        Loop_Batch_Txt;
    end
end


