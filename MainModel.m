%%% Script principal pour le lancer de batchs de tests
%%% Entrées : Chemins précisés ci-dessous (p et PathPreSet) menant aux données c3d et aux trajectoires cibles
%%% Sorties : Données issues de la boucle d'optimisation, sauvegardées en .mat dans SavePath
clear all;
close all;
clc;

flag = struct;

PathFiles = "";

model = struct;

% Load Model - Posture Description

model.description = load(strcat(PathFiles,'description.txt'));

% Load Model - Posture Reference

model.reference = load(strcat(PathFiles,'reference.txt'));

% Load Tracks - Empreintes Cibles 2 XYZ coordinates, in Pelvic frame coordinates 

model.prints = load(strcat(PathFiles,'prints.txt'));

% Load Initial Gait - Gait to be transformed / 11 * 100 vector of joint coordinates in splines 


% TODO Get Frames from Postures 



% TODO IK btwn Descr & Reference --> Ref State
% TODO IK btwn Ref & Initial posture (InitialGait(:,1)) --> Initial State
% TODO Inject variables into prep shit for LoopBatch













