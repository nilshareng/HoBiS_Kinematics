%%% Handles the difference between loading a .c3d walk and a simple model

model = struct;

p = 'C:\Users\nhareng\Documents\Prog\Matlab\Test\warping 2006\';
% Param = [Fem1g',Fem6g',Tal1g', Fem1d', Fem6d', Tal1d'];

% Load Model - Posture Description

model.description = load(strcat(p,'description.txt'));

% Load Model - Posture Reference

model.reference = load(strcat(p,'reference.txt'));
flag.steps = 0;

if flag.steps
    model.prints = load(strcat(p,'prints.txt'));
    X = model.prints;
end
% % Load Tracks - Empreintes Cibles 2 XYZ coordinates, in Pelvic frame coordinates 

% flag.prints = 0;
flag.prints = 1;


% Load Initial Gait - Gait to be transformed / 11 * 100 vector of joint coordinates in splines 

% PathPreSet = "";

% model.TA =  load(strcat(PathPreSet,"")) ;
InitialGaitPath = strcat(p,'poulaineAppui.txt');

if strcmp(InitialGaitPath(:,end-2:end),'mat')
    
SelectedPreset = 'antho012.mat';
load(strcat(PathPreSet,SelectedPreset));
Period = length(PN);

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

% Load Joint ranges
model.jointRangesMax =[-45; -45; -45; -90; -30; -60; 15; -90; -30; -60; 15]*pi/180;
model.jointRangesMin=[45; 45; 45; 20; 45; 30; 110; 20; 45; 30; 110]*pi/180;
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
% PosC = model.gait(1,:);
PosC = tmp;
deltaX = 0;
Angles = zeros(((max(size(model.gait)))+1)*10,11);
Angles(1,:) = anglesR2I';
NAngles = Angles(1,:);
d = 0.01;
thresh = 0.01;

%%
% Premiere IK pour positionner en Ref

NAngles = zeros(11,1);
Target= model.gait(1,:)';
Target = [Rmarkers.LTal1' ; Rmarkers.RTal1']/1000;
KinforMin(NAngles,Target,DParam,DReperes)
Angles = zeros(11,1);
options = optimset('Display','Iter','TolFun',1e-5);
AnglesDesc2Ref = fminsearch(@(Angles) KinforMin(Angles,Target,DParam,DReperes),NAngles,options);
CurrentPos = fcineshort(AnglesDesc2Ref,DParam,DReperes);
error = CurrentPos*1000 - [Rmarkers.LTal1' ; Rmarkers.RTal1'];%/1000
norm(error)
%%% Fonctionnel au dessus



% Deuxieme IK pour positionner en frame ini walk. LE reste fout la merde
TmpAngles = AnglesDesc2Ref;
thresh = 0.01;
options = optimset('TolFun',thresh);
mem = [];
steplength=0.01; %
globaltarget = model.gait(1,:)';
globalerror = [];
globalerror =  [globalerror , globaltarget - CurrentPos];

c = 0;
while c<100 && norm(globalerror(:,end))>thresh
c = c+1;
currenttarget = globaltarget-CurrentPos / norm(globaltarget-CurrentPos) * steplength;
delta = currenttarget - CurrentPos;
TmpAngles = pinv(JacShort(TmpAngles,DParam,DReperes)) * delta;
CurrentPos = fcineshort(TmpAngles,DParam, DReperes);
globalerror =  [globalerror , globaltarget - CurrentPos];
norm(globalerror(:,end))
end

%%%

for i = 1:1/steplength
% c=0;
% while norm(model.gait(1,:)' - CurrentPos)> thresh
%     c=c+1
    Target = (model.gait(1,:)' - CurrentPos);
    Target = CurrentPos + (Target*steplength / norm(Target))*i;
    TmpAngles = fminsearch(@(Angles) KinforMin(Angles,Target,DParam,DReperes),TmpAngles,options);
%     if KinforMin(TmpAngles,Target,DParam,DReperes)>thresh
%         Target = (CurrentPos + model.gait(1,:)')* i/50;
%         TmpAngles = fminsearch(@(Angles) KinforMin(Angles,Target,DParam,DReperes),TmpAngles,options);
%     end
    CurrentPos = fcineshort(TmpAngles,DParam,DReperes);
    mem = [mem, CurrentPos];
    if size(mem,2)>2 && norm(mem(:,end-1)-mem(:,end))<thresh/10
        warning('shit happened here')
        break;
    end
    error = [error, CurrentPos - Target];
    globalerror = [globalerror , (model.gait(1,:)' - CurrentPos)];
    norm(error(:,end))
end
AnglesRef2Ini1 = TmpAngles;
%%%
% Target = model.gait(1,:)';
% AnglesRef2Ini2 = fminsearch(@(Angles) KinforMin(Angles,Target,DParam,DReperes),AnglesRef2Ini1,options);
% CurrentPos = fcineshort(AnglesRef2Ini2,DParam,DReperes);
% CurrentPos - Target


%%%
test = zeros(11,1);
option=[1 200*10 0.03 2 0.5 0.002];
AnglesRef2Ini = MDS('KinforMin',test,AnglesDesc2Ref,Target,DParam,DReperes)
fcineshort(AnglesRef2Ini,DParam,DReperes) - model.gait(1,:)'


% IK sur le cycle de marche visé

x0 = AnglesDesc2Ref;
NewAngles = IK(model.Ngait,thresh,NAngles,RParam,RReperes,x0)';
figure;
plot(NewAngles);

NewPoul= [];
for i = 1:size(NewAngles,1)
    NewPoul =[NewPoul, fcineshort(NewAngles(i,:), NRParam, RReperes)];
end
NewPoul = NewPoul';
figure;
for i=1:6
    subplot(2,3,i)
    hold on;
    plot([NewPoul(:,i);NewPoul(1,i)]);
    plot(model.gait(:,i));
    plot(model.Ngait(:,i));
end

%%

for i = 1:(max(size(model.gait)))
    if i == (max(size(model.gait)))
        TargetX = model.gait(2,:)';
    else
        TargetX = model.gait(i+1,:)';
    end
%     deltaX = (TargetX - PosC);
%     delta = norm(deltaX)*d;
    
    deltaX = (TargetX - PosC);
    tmp = Angles(i,:);
    figure;
    hold on;
    while (norm(deltaX) > thresh)
        deltaX = (TargetX - PosC);
        delta = (deltaX/norm(deltaX)) * d;
        J = JacShort(NAngles, NRParam, RReperes);
        NAngles = pinv(J)*delta;
        PosC = fcineshort(NAngles, NRParam, RReperes);
        tmp = tmp + NAngles';
        plot(norm(deltaX));
    end
    i
%     for j = 1:10
%         deltaX = (TargetX - PosC);
%         deltax = deltaX/delta*d*j;%(1/100-j);
%         J = JacShort(NAngles, NRParam, RReperes);
%         NAngles = pinv(J) * deltax;
%         % TODO passer en fcineH3
%         %     PosC = fcine_numerique_H2(NAngles,Param(:,1)',Param(:,2)',Param(:,3)',Param(:,4)',Param(:,5)',Param(:,6)', R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local);
%         PosC = fcineshort(NAngles, NRParam, RReperes);
%         
%         Angles((i-1)*10+j,:) = NAngles';
%     end
    Angles(i+1,:) = tmp';
end
figure;
plot(Angles)
Angles = Angles(2:end-1,:);
rate = 60;

angles = [];
for i =10:10:600
    angles =[angles ;Angles(i,:)]; 
end
Angles = angles;
figure;
plot(Angles);
NewPoul = [];

for i = 1:size(Angles,1)
    NewPoul =[NewPoul, fcineshort(Angles(i,:), NRParam, RReperes)];
end
NewPoul = NewPoul';
figure;
for i=1:6
    subplot(2,3,i)
    hold on;
    plot([NewPoul(:,i);NewPoul(1,i)]);
%     plot(model.gait(:,i));
%     plot(model.Ngait(:,i));
end
% Approximation par spline des Traj Angulaires
PolA = [];
for i = 1:11
        [Var , PolT] = Curve2Spline(Angles(:,i));
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
    ComputedPoulaine = [ComputedPoulaine ; fcine_numerique_H2(NewCurve(i,:),Param(:,1)',Param(:,2)',Param(:,3)',Param(:,4)',Param(:,5)',Param(:,6)', R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local)'];
end

% NewSplines
PolP = [];
SplinedComputedPoulaine = [];
for i = 1:6
        [Var, PolT] = Curve2Spline(ComputedPoulaine(:,i));
        SplinedComputedPoulaine(:,i) = Var';
        PolP = [PolP ; [i*ones(size(PolT(:,2:end),1),1) , PolT(:,2:end)]];
end

% Affichages comparatif

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

GT = PN;
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