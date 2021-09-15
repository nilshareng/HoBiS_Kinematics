%%% Handles the difference between loading a .c3d walk and a simple model

model = struct;

PathFiles = "C:\Users\nhareng\Documents\Prog\Matlab\Test\warping 2006\";

% Load Model - Posture Description

model.description = load(strcat(PathFiles,'descriptionAstrid.txt'));

% Load Model - Posture Reference

model.reference = load(strcat(PathFiles,'referenceAstrid.txt'));

% Load Tracks - Empreintes Cibles 2 XYZ coordinates, in Pelvic frame coordinates 

model.prints = load(strcat(PathFiles,'prints.txt'));

% Load Initial Gait - Gait to be transformed / 11 * 100 vector of joint coordinates in splines 

PathPreSet = "";

model.TA =  load(strcat(PathPreSet,"")) ;
model.gait = load(strcat(PathFiles,"poulaineAppui.txt")) ;

% Load Joint ranges
model.jointRangesMax =[-45; -45; -45; -90; -30; -60; 15; -90; -30; -60; 15]*pi/180;
model.jointRangesMin=[45; 45; 45; 20; 45; 30; 110; 20; 45; 30; 110]*pi/180;

PathButees = "";

model.jointRanges = load(strcat(PathButees,""));
%%

% Compute Description to Reference angles 

anglesD2R = description2reference(strcat(PathFiles,'description.txt'),strcat(PathFiles,'reference.txt'));

% Compute Reference to initial trajectory angles

anglesR2I = position_initiale(anglesD2R,model.gait(1,1:3), model.gait(1,1:3),  model.jointRangesMax, model.jointRangesMin, Param(:,1)', Param(:,2)', Param(:,3)', Param(:,4)', Param(:,5)', Param(:,6))';


% 











%%%