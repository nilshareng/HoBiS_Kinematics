%%% Handles the difference between loading a .c3d walk and a simple model

% model = struct;

% p = 'C:\Users\nhareng\Documents\Prog\Matlab\Test\warping 2006\';
% Param = [Fem1g',Fem6g',Tal1g', Fem1d', Fem6d', Tal1d'];

% Load Model - Posture Description


% Load Model - Posture Reference

% % Load Tracks - Empreintes Cibles 2 XYZ coordinates, in Pelvic frame coordinates 

% flag.prints = 0;
flag.prints = 1;


% Load Initial Gait - Gait to be transformed / 11 * 100 vector of joint coordinates in splines 

% PathPreSet = "";

% model.TA =  load(strcat(PathPreSet,"")) ;
% InitialGaitPath = strcat(p,'poulaineAppui.txt');
MaxReachL = norm(Markers.LTal1 - Markers.LFem6) + norm(Markers.LFem6 - Markers.LHRC) + norm(Markers.LHRC - [0 0 0]);
MaxReachR = norm(Markers.RTal1 - Markers.RFem6) + norm(Markers.RFem6 - Markers.RHRC) + norm(Markers.RHRC - [0 0 0]);
MaxReach = (MaxReachL+MaxReachR)/2;

if strcmp(InitialGaitPath(:,end-2:end),'mat')
    
%     SelectedPreset = 'antho012.mat';
    tmpX = X;
    load(strcat(InitialGaitPath));
    PX = X;
    PX = Xtreat(PX);
    X = tmpX;
    [~, t] = sort(abs(min(PN(:,1:3))) + abs(max(PN(:,1:3))));
    t = [t , t+3];
    PN = PN(:,t);
    PN = [PN(:,1) -1*PN(:,2) PN(:,3) PN(:,4) -1*PN(:,5) PN(:,6)];
%     OPN = [OPN(:,1) -1*OPN(:,2) OPN(:,3) OPN(:,4) -1*OPN(:,5) OPN(:,6)];

%     [X(1:3,3:-1:1), PN] = Ratio2merde(PX,PN,X,OPN,MaxReach);    
%     X = Xtreat(X);
    
    InitialGait = PN;
    Period = max(size(InitialGait));
    midPeriod = fix((Period)/2);
    SplinedPoulaine = [];
    for i = 1:3
        SplinedPoulaine(:,i) = Curve2Spline(InitialGait(:,i));
    end
    % Symetrisation
    SplinedPoulaine = GaitSymetrisation(SplinedPoulaine);
    
    % Attribution
    model.gait = SplinedPoulaine;
    if flag.prints
       figure;
       hold on;
       title('Comparison loaded poulaine vs Spline poulaine');
       for i = 1:6
           subplot(2,3,i)
           hold on;
           plot(InitialGait(:,i));
           plot(SplinedPoulaine(:,i));
       end
    end

elseif strcmp(InitialGaitPath(:,end-2:end),'txt')

    InitialGait = load(strcat(InitialGaitPath));
    if size(InitialGait,2) == 7
       InitialGait = InitialGait(:,2:end); 
    end
    % Spline Approx of Initial Gait
    Period = max(size(InitialGait));
    midPeriod = fix((Period)/2);
    SplinedPoulaine = [];
    for i = 1:3
        SplinedPoulaine(:,i) = Curve2Spline(InitialGait(:,i));
    end
    % Symetrisation
    SplinedPoulaine = GaitSymetrisation(SplinedPoulaine);
    
    % Attribution
    model.gait = SplinedPoulaine;
    if flag.prints
       figure;
       hold on;
       title('Comparison loaded poulaine vs Spline poulaine');
       for i = 1:6
           subplot(2,3,i)
           hold on;
           plot(InitialGait(:,i));
           plot(SplinedPoulaine(:,i));
       end
    end
end

if exist('OPN')
    N=[];
    for i = 1:min(size(model.gait,1),size(OPN,1))
        N = [N ; norm(OPN(i,1:3)), norm(OPN(i,4:6))];
    end
    MaxPoul = max(max(N))*1000;
    PoulaineRatioOld = 1;
    if MaxPoul > MaxReach
        PoulaineRatioOld = MaxReach / (MaxPoul+100) ;
    end
    
    X(:,3:5) = X(:,3:5)*PoulaineRatioOld;
end

% Load Joint ranges
% model.jointRangesMax =[-45; -45; -45; -90; -30; -60; 15; -90; -30; -60; 15]*pi/180;
% model.jointRangesMin=[45; 45; 45; 20; 45; 30; 110; 20; 45; 30; 110]*pi/180;
Tmp = [];
for i = 1: size(model.gait,1)
    Tmp = [Tmp, [RReperes.Monde(1:3,1:3)^-1 * model.gait(i,1:3)' ; RReperes.Monde(1:3,1:3)^-1 * model.gait(i,4:6)']];
end
model.Ngait = Tmp';
% PathButees = "";
% model.jointRanges = load(strcat(PathButees,""));

% Compute Description to Reference angles 
% Test with NG function

% anglesD2R = description2reference(strcat(p,'description.txt'),strcat(p,'reference.txt'));
% PosCheckRef = fcine_numerique_H2(anglesD2R,Param(:,1)',Param(:,2)',Param(:,3)',Param(:,4)',Param(:,5)',Param(:,6)', R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local);

% Compute Reference to initial trajectory angles
% Test with NG function

% anglesR2I = position_initiale(anglesD2R,model.gait(1,1:3), model.gait(1,1:3),  model.jointRangesMax, model.jointRangesMin, Param(1,:), Param(2,:), Param(:,3)', Param(4,:), Param(5,:), Param(6,:))';
% Vérif TA 
% PosCheckIni = fcine_numerique_H2(anglesR2I,Param(:,1)',Param(:,2)',Param(:,3)',Param(:,4)',Param(:,5)',Param(:,6)', R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local);
% PosCheckIni = fcineshort(anglesR2I, NRParam, RReperes)
anglesR2I = zeros(11,1);

% Compute full Initial Articular Trajectories corresponding to Initial Gait - Full IK, puis Spline 
% IK :

J = [];
PosC = model.gait(1,:);
% PosC = a;
deltaX = 0;
Angles = zeros(((max(size(model.gait)))+1)*10,11);
Angles(1,:) = anglesR2I';
NAngles = Angles(1,:);
d = 0.01;
thresh = 0.01;

%%
% Premiere IK pour positionner en Ref
model.gait = SplinedPoulaine;
NAngles = zeros(11,1);
Target= model.gait(1,:)';
Target = [Rmarkers.LTal1' ; Rmarkers.RTal1']/1000;
Angles = zeros(11,1);
KinforMin(Angles,Sequence,Target,Markers,Reperes)
options = optimset('Display','Iter','TolFun',1e-5);
AnglesDesc2Ref = fminsearch(@(Angles) KinforMin(Angles,Sequence,Target,Markers,Reperes),NAngles,options);
[CurrentPos, Cmarkers, Creperes] = fcinematique(AnglesDesc2Ref,Sequence,Markers,RReperes);
DisplayMarkers(Cmarkers,1,Creperes);
error = CurrentPos*1000 - [Rmarkers.LTal1' ; Rmarkers.RTal1'];%/1000
norm(error)
%%% Fonctionnel au dessus



% Deuxieme IK pour positionner en frame ini walk. LE reste fout la merde
TmpAngles = zeros(1,11);
thresh = 0.01;
options = optimset('TolFun',thresh);
mem = [];
steplength=0.01; %
globaltarget = [model.gait(1,4:6) , model.gait(1,1:3)]';
% globaltarget = [model.gait(1,1:3) , model.gait(1,4:6)]';
globalerror = [];
globalerror =  [globalerror , globaltarget - CurrentPos];
[CurrentPos, Cmarkers, Creperes] = fcinematique(TmpAngles,Sequence,Markers,RReperes);

TmpMarkers = Markers;
TmpMarkers.LTarget = globaltarget(1:3)*1000;
TmpMarkers.RTarget = globaltarget(4:6)*1000;
currenttarget = CurrentPos + (globaltarget-CurrentPos) / norm(globaltarget-CurrentPos) * steplength;
TmpMarkers.CRTarget = currenttarget(4:6)*1000;
TmpMarkers.CLTarget = currenttarget(1:3)*1000;
c = 0;

% while c<100 && norm(globalerror(:,end))>thresh
% c = c+1;
% currenttarget = CurrentPos + (globaltarget-CurrentPos) / norm(globaltarget-CurrentPos) * steplength;
% deltaX = (currenttarget - CurrentPos);
% deltatheta =  pinv(Jcinematique(TmpAngles,Sequence,TmpMarkers,Reperes)) * deltaX;
% TmpAngles = TmpAngles + deltatheta;
% [CurrentPos, TmpMarkers, TmpReperes] = fcinematique(TmpAngles,Sequence,TmpMarkers,Reperes);
% globalerror =  [globalerror , globaltarget - CurrentPos];
% norm(globalerror(:,end))
% TmpMarkers.LTarget = globaltarget(1:3)*1000;
% TmpMarkers.RTarget = globaltarget(4:6)*1000;
% TmpMarkers.CRTarget = currenttarget(1:3)*1000;
% TmpMarkers.CLTarget = currenttarget(4:6)*1000;
% DisplayMarkers(TmpMarkers,0,TmpReperes);
% view([45 -45 0])
% pause;
% close all;
% end
% DisplayMarkers(TmpMarkers,0,TmpReperes);

options = optimset('TolFun',1e-5);
NewAngles = [];
NewPoul = [];
% Markers.LPoul
% Markers.RPoul = 


N=[]; 
for i = 1:size(model.gait,1)
    N = [N ; norm(model.gait(i,1:3)), norm(model.gait(i,4:6))];
end
MaxPoul = max(max(N))*1000;
PoulaineRatio = 1;
if MaxPoul > MaxReach
    PoulaineRatio = MaxReach / (MaxPoul+100) ;
end

% close all;

% PoulaineRatio = 0.8;
% PoulaineRatio = 1;

model.gait = model.gait * PoulaineRatio;


% model.gait(:,3:3:6) = model.gait(:,3:3:6) + ones(size(model.gait(:,3:3:6))) * (-0.2);
NewPoul=[];

GaitMarkers = [];
% GaitMarkers = struct;
s = 10;
%%% "IK" Disgusting fminseach instead, cause it (kinda) works 
for i =1:size(model.gait,1)
    globaltarget = [model.gait(i,1:3) , model.gait(i,4:6)]';
%     tmp = fminsearch(@(Angles) KinforMin(Angles,Sequence,globaltarget,Markers,Reperes),TmpAngles,options);
    tmp = TmpAngles;
    for j = 1:s
%         localtarget = CurrentPos + ((globaltarget - CurrentPos)/norm(globaltarget - CurrentPos)) /10;
        localtarget = CurrentPos + ((globaltarget - CurrentPos)/(s+1-j));
        deltaX = CurrentPos - localtarget;
        J = Jcinematique(TmpAngles,Sequence,Markers,RReperes);
        Jp = pinv(J);
        deltatheta =  Jp * deltaX; 
        proj = Jp*J;
        tmp = tmp + deltatheta';
        deltatheta2 = ones(1,11)*0.001;
        delta2 = fminsearch(@(deltatheta2) ArticularCost(tmp,deltatheta2,proj,model.jointRangesMin, model.jointRangesMax),deltatheta2);
        tmp = tmp + (proj*delta2')';
        CurrentPos = fcinematique(tmp,Sequence,Markers,RReperes);
    end
    NewAngles = [NewAngles ; tmp];
    [CurrentPos, Cmarkers, Creperes] = fcinematique(tmp,Sequence,Markers,RReperes);
    Cmarkers.CLTarget = globaltarget(1:3)'*1000;
    Cmarkers.CRTarget = globaltarget(4:6)'*1000;
    Cmarkers.LPoul = model.gait(:,1:3)*1000;
    Cmarkers.RPoul = model.gait(:,4:6)*1000;
    
    GaitMarkers = [GaitMarkers,Cmarkers];
    
    NewPoul = [NewPoul ; CurrentPos'];
%     DisplayMarkers(Cmarkers,1,Creperes);
%     pause(1/60); 
    
    
%     pause;
%     close all;
end
for i =1:size(model.gait,1)
    GaitMarkers(i).LOPoul = NewPoul(:,1:3)*1000;
    GaitMarkers(i).ROPoul = NewPoul(:,4:6)*1000;
    GaitMarkers(i).X = X(:,3:5)*1000;
end
DisplayCurves(NewPoul);
% figure;
% hold on;
DisplayCurves(NewAngles);

%%% Filtering the noisy results
freq=5;
[b,a] = butter(2 , freq/(0.5*Period) , 'low');
S = size(NewAngles, 1);
NewAnglesF=filtfilt(b,a,[NewAngles ; NewAngles ; NewAngles]);
GaitMarkers2 = [];
NewPoul2 = [];

for i = 1:Period
    [tmp, tmpMarkers, tmpReperes] = fcinematique(NewAnglesF(S+i,:),Sequence,Markers,RReperes);
    GaitMarkers2 = [GaitMarkers2,tmpMarkers];
    
    NewPoul2 = [NewPoul2 ; tmp'];
end
    
[b,a] = butter(2 , freq/(0.5*Period) , 'low');
NewPoulF=filtfilt(b,a,[NewPoul; NewPoul; NewPoul]);
rate = Period;
%%% Display 
f=figure;
hold on;
model.invgait = [model.gait(:,4:6) , model.gait(:, 1:3)];
for i = 1:6
    subplot(2,3,i)
    hold on;
    plot(NewPoulF(size(NewPoul,1)+1:2*size(NewPoul,1),i));
    plot(NewPoul2(:,i));
    plot(model.gait(:,i));
end

figure;
hold on;
for i = 1:11
    subplot(4,3,i)
    plot(NewAnglesF(S+1:2*S,i));
end

NewPoulF = NewPoulF(size(NewPoul,1)+1:2*size(NewPoul,1),:);
% DisplayGait(GaitMarkers);

NewAnglesF = NewAnglesF(S+1:2*S,:);

% for i = 1:size(GaitMarkers.RFWT,1)
%     
%     DisplayMarkers(,2)
% end
%%%
%%
% %%
% for i = 1:1/steplength
% % c=0;
% % while norm(model.gait(1,:)' - CurrentPos)> thresh
% %     c=c+1
%     Target = (model.gait(1,:)' - CurrentPos);
%     Target = CurrentPos + (Target*steplength / norm(Target))*i;
%     TmpAngles = fminsearch(@(Angles) KinforMin(Angles,Target,DParam,DReperes),TmpAngles,options);
% %     if KinforMin(TmpAngles,Target,DParam,DReperes)>thresh
% %         Target = (CurrentPos + model.gait(1,:)')* i/50;
% %         TmpAngles = fminsearch(@(Angles) KinforMin(Angles,Target,DParam,DReperes),TmpAngles,options);
% %     end
%     CurrentPos = fcineshort(TmpAngles,DParam,DReperes);
%     mem = [mem, CurrentPos];
%     if size(mem,2)>2 && norm(mem(:,end-1)-mem(:,end))<thresh/10
%         warning('shit happened here')
%         break;
%     end
%     error = [error, CurrentPos - Target];
%     globalerror = [globalerror , (model.gait(1,:)' - CurrentPos)];
%     norm(error(:,end))
% end
% AnglesRef2Ini1 = TmpAngles;
% %%%
% % Target = model.gait(1,:)';
% % AnglesRef2Ini2 = fminsearch(@(Angles) KinforMin(Angles,Target,DParam,DReperes),AnglesRef2Ini1,options);
% % CurrentPos = fcineshort(AnglesRef2Ini2,DParam,DReperes);
% % CurrentPos - Target
% 
% 
% %%%
% test = zeros(11,1);
% option=[1 200*10 0.03 2 0.5 0.002];
% AnglesRef2Ini = MDS('KinforMin',test,AnglesDesc2Ref,Target,DParam,DReperes)
% fcineshort(AnglesRef2Ini,DParam,DReperes) - model.gait(1,:)'
% 
% 
% % IK sur le cycle de marche visé
% 
% x0 = AnglesDesc2Ref;
% NewAngles = IK(model.Ngait,thresh,NAngles,RParam,RReperes,x0)';
% figure;
% plot(NewAngles);
% 
% NewPoul= [];
% for i = 1:size(NewAngles,1)
%     NewPoul =[NewPoul, fcineshort(NewAngles(i,:), NRParam, RReperes)];
% end
% NewPoul = NewPoul';
% figure;
% for i=1:6
%     subplot(2,3,i)
%     hold on;
%     plot([NewPoul(:,i);NewPoul(1,i)]);
%     plot(model.gait(:,i));
%     plot(model.Ngait(:,i));
% end
% 
% %
% 
% for i = 1:(max(size(model.gait)))
%     if i == (max(size(model.gait)))
%         TargetX = model.gait(2,:)';
%     else
%         TargetX = model.gait(i+1,:)';
%     end
% %     deltaX = (TargetX - PosC);
% %     delta = norm(deltaX)*d;
%     
%     deltaX = (TargetX - PosC);
%     tmp = Angles(i,:);
%     figure;
%     hold on;
%     while (norm(deltaX) > thresh)
%         deltaX = (TargetX - PosC);
%         delta = (deltaX/norm(deltaX)) * d;
%         J = JacShort(NAngles, NRParam, RReperes);
%         NAngles = pinv(J)*delta;
%         PosC = fcineshort(NAngles, NRParam, RReperes);
%         tmp = tmp + NAngles';
%         plot(norm(deltaX));
%     end
%     i
% %     for j = 1:10
% %         deltaX = (TargetX - PosC);
% %         deltax = deltaX/delta*d*j;%(1/100-j);
% %         J = JacShort(NAngles, NRParam, RReperes);
% %         NAngles = pinv(J) * deltax;
% %         % TODO passer en fcineH3
% %         %     PosC = fcine_numerique_H2(NAngles,Param(:,1)',Param(:,2)',Param(:,3)',Param(:,4)',Param(:,5)',Param(:,6)', R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local);
% %         PosC = fcineshort(NAngles, NRParam, RReperes);
% %         
% %         Angles((i-1)*10+j,:) = NAngles';
% %     end
%     Angles(i+1,:) = tmp';
% end
% figure;
% plot(Angles)
% Angles = Angles(2:end-1,:);
% rate = 60;
% 
% angles = [];
% for i =10:10:600
%     angles =[angles ;Angles(i,:)]; 
% end
% Angles = angles;
% figure;
% plot(Angles);
% NewPoul = [];
% 
% for i = 1:size(Angles,1)
%     NewPoul =[NewPoul, fcineshort(Angles(i,:), NRParam, RReperes)];
% end
%%
close all;
MarkersReference = Markers;
% NewPoul = NewPoul';
figure;
for i=1:6
    subplot(2,3,i)
    hold on;
    plot([NewPoulF(:,i);NewPoulF(1,i)]);
%     plot(model.gait(:,i));
%     plot(model.Ngait(:,i));
end
% Approximation par spline des Traj Angulaires
PolA = [];
SplinedAngles = [];
for i = 1:11
        [Var , PolT] = Curve2Spline(NewAnglesF(:,i));
        SplinedAngles(:,i) = Var';
        PolA = [PolA ; [i*ones(size(PolT(:,2:end),1),1) , PolT(:,2:end)]];
end
% Symetrisation
SplinedAngles = GaitSymetrisation(SplinedAngles);

% PCA

PCA = ForgePCA(SplinedAngles,0:1/(rate-1):1 ,1 );
PolA = [];
for i =1:PCA(end,1)
    PolA = [PolA ; PC_to_spline(PCA(PCA(:,1)==i,2:3),1)];
end
PolA(:,1) = PCA(:,1);

NewCurve = [];
% for i =1:Pol(end,1)
%     NewCurve = [NewCurve , spline_to_curve_int(Pol(Pol(:,1)==i,:),Period)'];
% end
for i =1:Period+1
    tmp=zeros(1,11);
    for j=1:11
        tmp(j) = EvalSpline(PolA(PolA(:,1)==j,:),((i-1))/Period);
    end
    NewCurve = [NewCurve ; tmp];
end
% Expression de la poulaine correspondante à l'approximation TA
ComputedPoulaine =[];
for i = 1:max(size(SplinedAngles))
    ComputedPoulaine = [ComputedPoulaine ; fcinematique(SplinedAngles(i,:),Sequence, MarkersReference, RReperes)'];
end

% NewSplines
PolP = [];
SplinedComputedPoulaine = [];
for i = 1:6
        [Var, PolT] = Curve2Spline(ComputedPoulaine(:,i));
        SplinedComputedPoulaine(:,i) = Var';
        PolP = [PolP ; [i*ones(size(PolT(:,2:end),1),1) , PolT(:,2:end)]];
end

tmpPCA = ForgePCA(SplinedAngles,0:1/(rate-1):1 ,1 );
Cheat = tmpPCA(tmpPCA(:,1)<=2,:);
Cheat(2:2:end,2) = Cheat(1:2:end,2) + 0.5;
Cheat(2:2:end,3) = -1*Cheat(1:2:end,3);
tmpPCA(tmpPCA(:,1)<=2,:) = Cheat;
NPCA = [tmpPCA(tmpPCA(:,1)<=2,:) ; PCA(PCA(:,1)>2,:)];

NPolA = [];
for i = 1:11
    temp = PC_to_spline(NPCA(NPCA(:,1)==i,2:3),1);
    NPolA = [NPolA ; i*ones(size(temp,1),1), temp(:,2:end)];
end

[tmpP2, tmpTA2] = Sampling_txt(NPolA,Period,Sequence,Markers,Reperes);

figure(6);
hold on;
DisplayCurves(tmpTA2,6);
DisplayCurves(SplinedAngles,6);

figure(2);
hold on;
DisplayCurves(tmpP2,2);
DisplayCurves(ComputedPoulaine,2);
for i = 1   :3
    subplot(3,3,X(i,1));
    hold on;
    plot(X(i,2)*(size(PN,1)-2)+1,X(i,X(i,1)+2),'kx');
end
for i = 4:6
    subplot(3,3,X(i,1));
    hold on;
    plot(X(i,2)*(size(PN,1)-2)+1,X(i,X(i,1)-1),'kx');
end








%% Affichages comparatif

if flag.prints
       figure;
       hold on;
       title('Comparison loaded poulaine vs Splined poulaine');
       for i = 1:6
           subplot(2,3,i)
           hold on;
           plot(InitialGait(:,i));
           plot(SplinedPoulaine(:,i));
           plot(ComputedPoulaine(:,i));
           plot(SplinedComputedPoulaine(:,i));
       end
end


% Setting up variables to plug into the PreLoop / Loop
PN = SplinedComputedPoulaine;

GT = model.gait;
Pol = PolA;
mid = fix(max(size(SplinedComputedPoulaine,1))/2);
NewCurve = SplinedAngles;
% for i = 1:max(size(SplinedAngles),1)
%     PCA = [PolA(i,1:2), EvalSpline(PolA(i,:),PolA(i,2)) ]; % PC : N X T : idNumber, Xvalue, Time value
% end
% for i = 1:max(size(SplinedComputedPoulaine))
%     PCP = [PolP(i,1:2), ]; % PC : N X T : idNumber, Xvalue, Time value
% end
% 





















%%%