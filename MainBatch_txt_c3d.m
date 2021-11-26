%%% M�J 06/09/2021
%%% Script principal pour le lancer de batchs de tests
%%% Entr�es : Chemins pr�cis�s ci-dessous (p et PathPreSet) menant aux donn�es c3d et aux trajectoires cibles
%%% Sorties : Donn�es issues de la boucle d'optimisation, sauvegard�es en .mat dans SavePath
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

definput = {'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\description.txt',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\reference.txt',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\NewestPresets\',...
    '',...
    '70',...
    '[45 45 45 20 45 30 110 20 45 30 110]',...
    '[-45 -45 -45 -90 -30 -60 15 -90 -30 -60 15]',...
    'C:\Users\nhareng\Desktop\CodeCommente\hobis\Resultats\Txt'};


RPath = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\NewestPresets\';
Ressources = ls(RPath);
% mkdir(strcat(definput{end},'\Batch1'));

counter1 = 1;
counter2 = 0;
FileSelect1 = [1 1 1 3 ...% Poulaines
    12 ...
    12 12 12 14 ...
    23 ...
    23 23 23 25 ...
    34 ...
    ] + 3;
FileSelect2 = [1 3 12 1 ... %Empreintes
    1 ...
    12 14 23 12 ...
    12 ... 
    23 25 34 23 ...
    23 ...
    ] + 3;
FileNames = {...
    'P_antho012_E_antho012', 'P_antho012_E_antho028' , 'P_antho012_E_armel012' , 'P_antho028_E_antho012', ...
    'P_armel012_E_antho012',...
    'P_armel012_E_armel012', 'P_armel012_E_armel028', 'P_armel012_E_hassane012', 'P_armel028_E_armel012', ...
    'P_hassane012_E_armel012',...
    'P_hassane012_E_hassane012', 'P_hassane012_E_hassane028' , 'P_hassane012_E_laurent012', 'P_hassane028_E_hassane012'...
    'P_laurent012_E_hassane012', ...
    };

for iii = 2:size(FileSelect1,2)
    % Empreintes loading
    load(strcat(RPath, Ressources(FileSelect2(iii),:)));
    [~, t] = sort(abs(min(PN(:,1:3))) + abs(max(PN(:,1:3))));
    t = [t , t+3];
    PN = PN(:,t);
    PN = [PN(:,1) -1*PN(:,2) PN(:,3) PN(:,4) -1*PN(:,5) PN(:,6)];
    OPN = PN;
    if FileSelect1(iii)==FileSelect2(iii)
        X(1:3,3:5) = X(1:3,3:5) + rand(size(X(1:3,3:5))) * 0.01;
    end
    X =Xtreat(X);
    close all;
    clc;
    clearvars -except definput flag p RPath Ressources iii X OX OPN FileNames FileSelect1 FileSelect2
    % Poulaine loading
    definput{5} = strcat(definput{4},Ressources(FileSelect1(iii),:));
    Automatized_Batch_Scripted;
    save(strcat(SavePath,'\Batch3\',FileNames{iii},'Results.mat'),'Results');
end

