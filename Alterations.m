%%% Handles the difference between loading a .c3d walk and a simple model

model = struct;

PathFiles = 'C:\Users\nhareng\Documents\Prog\Matlab\Test\warping 2006\';
Param = [Fem1g',Fem6g',Tal1g', Fem1d', Fem6d', Tal1d'];

% Load Model - Posture Description

model.description = load(strcat(PathFiles,'descriptionAstrid.txt'));

% Load Model - Posture Reference

model.reference = load(strcat(PathFiles,'referenceAstrid.txt'));

% % Load Tracks - Empreintes Cibles 2 XYZ coordinates, in Pelvic frame coordinates 

% model.prints = load(strcat(PathFiles,'prints.txt'));

% Load Initial Gait - Gait to be transformed / 11 * 100 vector of joint coordinates in splines 

% PathPreSet = "";

% model.TA =  load(strcat(PathPreSet,"")) ;
InitialGaitPath = strcat(PathFiles,'poulaineAppui.txt');

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

% PathButees = "";
% model.jointRanges = load(strcat(PathButees,""));

% Compute Description to Reference angles 
% Test with NG function

anglesD2R = description2reference(strcat(PathFiles,'description.txt'),strcat(PathFiles,'reference.txt'));
PosCheckRef = fcine_numerique_H2(anglesD2R,Param(:,1)',Param(:,2)',Param(:,3)',Param(:,4)',Param(:,5)',Param(:,6)', R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local);

% Compute Reference to initial trajectory angles
% Test with NG function

anglesR2I = position_initiale(anglesD2R,model.gait(1,1:3), model.gait(1,1:3),  model.jointRangesMax, model.jointRangesMin, Param(:,1)', Param(:,2)', Param(:,3)', Param(:,4)', Param(:,5)', Param(:,6))';
% Vérif TA 
PosCheckIni = fcine_numerique_H2(anglesR2I,Param(:,1)',Param(:,2)',Param(:,3)',Param(:,4)',Param(:,5)',Param(:,6)', R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local);

% Compute full Initial Articular Trajectories corresponding to Initial Gait - Full IK, puis Spline 
% IK :
J = [];
% PosC = model.gait(1,:);
PosC = AnglesCheckIni;
deltaX = 0;
Angles = anglesR2I;

for i = 1:(max(size(model.gait))+1)
    
    if i == (max(size(model.gait))+1)
        TargetX = model.gait(1,:);
    else
        TargetX = model.gait(i+1,:)';
    end
    
    deltaX = TargetX - PosC;
    
    J = Jac_H2(Angles(end,:),Param(:,1)',Param(:,2)',Param(:,3)',Param(:,4)',Param(:,5)',Param(:,6)', R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local);
    
    NAngles = pinv(J) * deltaX;
    
    PosC = fcine_numerique_H2(NAngles,Param(:,1)',Param(:,2)',Param(:,3)',Param(:,4)',Param(:,5)',Param(:,6)', R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local);
    Angles = [Angles ; NAngles];
end

rate = 60;

% Approximation par spline des Traj Angulaires

for i = 1:11
        SplinedAngles(:,i) = Curve2Spline(Angles(:,i));
end
% Symetrisation
SplinedAngles = GaitSymetrisation(SplinedAngles);

% Expression de la poulaine correspondante à l'approximation TA
ComputedPoulaine
for i = 1:max(size(SplinedAngles))
    ComputedPoulaine = fcine
end
% 






















%%%