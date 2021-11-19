%%% MàJ 06/09/2021
%%% Script principal pour le lancer de batchs de tests
%%% Entrées : Chemins précisés ci-dessous (p et PathPreSet) menant aux données c3d et aux trajectoires cibles
%%% Sorties : Données issues de la boucle d'optimisation, sauvegardées en .mat dans SavePath
clear all;
close all;
clc;

p = genpath('.');
addpath(p);
%addpath("C:\Users\nhareng\Desktop\CodeCommente\hobis");

flag = struct;

flag.c3d = 0;
flag.txt = 0;
flag.dyn = 0;
flag.prints =0;

%% Input part
%Antho012 - Antho056

% prompt = ...
%     {'Enter Ressources directory (absolute) path:',...
%     'Enter Description Posture file (absolute) path:',...
%     'Enter Reference Posture file (absolute) path:',...
%     'Enter PreSets files directory (absolute) path:',...
%     'Initial Poulaine (absolute) path:',...
%     'Model Mass:', ...
%     'Model Articular Range of Motion Min Boundaries in degrees (Pelv*3, LHip*3, LKnee*1, RHip*3, RKnee*1):',...
%     'Model Articular Range of Motion Max Boundaries in degrees (Pelv*3, LHip*3, LKnee*1, RHip*3, RKnee*1):',...
%     'Footprints (target) left footstrike (gait% X Y Z):',...
%     'Footprints (target) left toe off (gait% X Y Z):',...
%     'Footprints (target) left min Z (gait% X Y Z):'...
%     'Enter Save folder Path:'};
% %     'Footprints (target) right footstrike (gait% X Y Z):',...
% %     'Footprints (target) right toe off (gait% X Y Z):',...
% %     'Footprints (target) right min Z (gait% X Y Z):'...
% dlgtitle = 'Inputs';
% dims = [1 110];

definput = {'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\description.txt',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\reference.txt',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\NewestPresets\',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\NewestPresets\antho012.mat',...
    '70',...
    '[45 45 45 20 45 30 110 20 45 30 110]',...
    '[-45 -45 -45 -90 -30 -60 15 -90 -30 -60 15]',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Resultats\Txt'};


RPath = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\NewestPresets\';
Ressources = ls(RPath);
mkdir(strcat(definput{end},'\Batch1'));
for iii = 6:11:(size(Ressources,1))
    for jjj = 4:(size(Ressources,1))
        if jjj~=iii
            load(strcat(RPath, Ressources(iii,:)));
            X = [X(1:3,2:5) ; X(1:3,2:5)];
            X = [X(:,1) , X(:,3) , -1*X(:,2) , X(:,4)];
            X(4:6,:) = [X(4:6,1) + 0.5 , X(4:6,2)*-1 , X(4:6,3:4)];
            X(X(1:6,1)>1,1) = X(X(1:6,1)>1,1) -1;
            X = [[2;2;3;5;5;6] , X];
            OPN = PN;
            close all;
            clc;
            clearvars -except definput flag p RPath Ressources iii jjj X OPN
            definput{5} = strcat(RPath,Ressources(iii,:));
            MainBatch_txt_c3d_Scripted;
            save(strcat(SavePath,'\Batch1\',num2str(fix(clock)),'.mat'),'Results');
            savefig(fConv,strcat(SavePath,'\Batch1\',num2str(fix(clock)),'.fig'));
            savefig(fFGait,strcat(SavePath,'\Batch1\',num2str(fix(clock)),'.fig'));
        end
    end
end


