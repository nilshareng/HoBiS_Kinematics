function [Distances,Stats] = PostTreatment(Path,Display)
% Stats and displays of results from previous algo
close all;
Ressources = ls(Path);
Ressources = Ressources(5:end,:); % Sale, valable pour \Batch3\
p = genpath('C:\Users\nhareng\Desktop\CodeCommente\hobis\');
addpath(p);
S = size(Ressources,1);

% PoulainesCollector
% TACollector

mkdir(Path,'Figures')

i = 0;
% while i < S-2
%     i = i+1;
%     C = rem(i,5);
%     switch C
%         case 1 % 12/12
%             load(strcat(Path, Ressources(i,:)));
%             D1 = Results;
%             
%             D1.PoulaineEmpreintesInput(:,2) = -1* D1.PoulaineEmpreintesInput(:,2);
%             D1.PoulaineEmpreintesInput(:,5) = -1* D1.PoulaineEmpreintesInput(:,5);
%             D1.PoulaineEmpreintesInputScaled(:,2) = -1* D1.PoulaineEmpreintesInputScaled(:,2);
%             D1.PoulaineEmpreintesInputScaled(:,5) = -1* D1.PoulaineEmpreintesInputScaled(:,5);
%             
%             DisplayCurves(D1.InitialPoulaine,i,'--');
%             DisplayCurves(D1.PoulaineEmpreintesInput,i,'--');
%             DisplayCurves(D1.PoulaineEmpreintesInputScaled,i,'--');
%             DisplayCurves(D1.FinalPoulaine,i);
%             
%             sgtitle(Ressources(i,1:end-4));
%             
%             f = figure(i);
%             hold on;
%             legend('Poulaine Initiale Scalée','Poulaine Empreintes','Poulaine Empreintes Scalée','Poulaine Finale');
%             DisplayX(D1.OriginalX, size(D1.PoulaineEmpreintesInput,1),i,'go');
%             DisplayX(D1.EmpreintesScaled, size(D1.PoulaineEmpreintesInput,1),i,'gx');
%             
%             subplot(3,3,7);
%             plot(D1.Convergence(1,:));
%             subplot(3,3,8);
%             plot(D1.Convergence(2,:));
% %             DisplayGait(Results.GaitMarkers,100);
%             
%         case 2 % 12/28 - 28/12
%             load(strcat(Path, Ressources(i,:)));
%             D1 = Results;
%             
%             D1.PoulaineEmpreintesInput(:,2) = -1* D1.PoulaineEmpreintesInput(:,2);
%             D1.PoulaineEmpreintesInput(:,5) = -1* D1.PoulaineEmpreintesInput(:,5);
%             D1.PoulaineEmpreintesInputScaled(:,2) = -1* D1.PoulaineEmpreintesInputScaled(:,2);
%             D1.PoulaineEmpreintesInputScaled(:,5) = -1* D1.PoulaineEmpreintesInputScaled(:,5);
%             
%             load(strcat(Path, Ressources(i+2,:)));
%             D2 = Results;
%             
%             D2.PoulaineEmpreintesInput(:,2) = -1* D2.PoulaineEmpreintesInput(:,2);
%             D2.PoulaineEmpreintesInput(:,5) = -1* D2.PoulaineEmpreintesInput(:,5);
%             D2.PoulaineEmpreintesInputScaled(:,2) = -1* D2.PoulaineEmpreintesInputScaled(:,2);
%             D2.PoulaineEmpreintesInputScaled(:,5) = -1* D2.PoulaineEmpreintesInputScaled(:,5);
%             
%             DisplayCurves(D1.InitialPoulaine,i,'--');
%             DisplayCurves(D1.PoulaineEmpreintesInput,i,'--');
%             DisplayCurves(D1.PoulaineEmpreintesInputScaled,i,'--');
%             DisplayCurves(D1.FinalPoulaine,i);
%             
%             DisplayCurves(D2.InitialPoulaine,i,'--');
%             DisplayCurves(D2.PoulaineEmpreintesInput,i,'--');
%             DisplayCurves(D2.PoulaineEmpreintesInputScaled,i,'--');
%             DisplayCurves(D2.FinalPoulaine,i);
%             
%             sgtitle(Ressources(i,1:end-4));
%             legend('Poulaine Initiale Scalée1','Poulaine Empreintes1','Poulaine Empreintes Scalée1','Poulaine Finale1','Poulaine Initiale Scalée2','Poulaine Empreintes2','Poulaine Empreintes Scalée2','Poulaine Finale2');
%             
%             f = figure(i);
%             hold on;
%             DisplayX(D1.OriginalX, size(D1.PoulaineEmpreintesInput,1),i,'go');
%             DisplayX(D1.EmpreintesScaled, size(D1.PoulaineEmpreintesInput,1),i,'gx');
%             DisplayX(D2.OriginalX, size(D2.PoulaineEmpreintesInput,1),i,'ro');
%             DisplayX(D2.EmpreintesScaled, size(D2.PoulaineEmpreintesInput,1),i,'rx');
%             
%             subplot(3,3,7);
%             hold on;
%             plot(D1.Convergence(1,:));
%             plot(D2.Convergence(1,:));
%             subplot(3,3,8);
%             hold on;
%             plot(D1.Convergence(2,:));
%             plot(D2.Convergence(2,:));
%             
%         case 3 % 12/12 - 12/12
%             load(strcat(Path, Ressources(i,:)));
%             D1 = Results;
%             
%             D1.PoulaineEmpreintesInput(:,2) = -1* D1.PoulaineEmpreintesInput(:,2);
%             D1.PoulaineEmpreintesInput(:,5) = -1* D1.PoulaineEmpreintesInput(:,5);
%             D1.PoulaineEmpreintesInputScaled(:,2) = -1* D1.PoulaineEmpreintesInputScaled(:,2);
%             D1.PoulaineEmpreintesInputScaled(:,5) = -1* D1.PoulaineEmpreintesInputScaled(:,5);
%             
%             load(strcat(Path, Ressources(i+2,:)));
%             D2 = Results;
%             
%             D2.PoulaineEmpreintesInput(:,2) = -1* D2.PoulaineEmpreintesInput(:,2);
%             D2.PoulaineEmpreintesInput(:,5) = -1* D2.PoulaineEmpreintesInput(:,5);
%             D2.PoulaineEmpreintesInputScaled(:,2) = -1* D2.PoulaineEmpreintesInputScaled(:,2);
%             D2.PoulaineEmpreintesInputScaled(:,5) = -1* D2.PoulaineEmpreintesInputScaled(:,5);
%             
%             DisplayCurves(D1.InitialPoulaine,i,'--');
%             DisplayCurves(D1.PoulaineEmpreintesInput,i,'--');
%             DisplayCurves(D1.PoulaineEmpreintesInputScaled,i,'--');
%             DisplayCurves(D1.FinalPoulaine,i);
%                         
%             DisplayCurves(D2.InitialPoulaine,i,'--');
%             DisplayCurves(D2.PoulaineEmpreintesInput,i,'--');
%             DisplayCurves(D2.PoulaineEmpreintesInputScaled,i,'--');
%             DisplayCurves(D2.FinalPoulaine,i);
%             
%             sgtitle(Ressources(i,1:end-5));
%             legend('Poulaine Initiale Scalée1','Poulaine Empreintes1','Poulaine Empreintes Scalée1','Poulaine Finale1','Poulaine Initiale Scalée2','Poulaine Empreintes2','Poulaine Empreintes Scalée2','Poulaine Finale2');
%             
%             f = figure(i);
%             hold on;
%             DisplayX(D1.OriginalX, size(D1.PoulaineEmpreintesInput,1),i,'go');
%             DisplayX(D1.EmpreintesScaled, size(D1.PoulaineEmpreintesInput,1),i,'gx');
%             DisplayX(D2.OriginalX, size(D2.PoulaineEmpreintesInput,1),i,'ro');
%             DisplayX(D2.EmpreintesScaled, size(D2.PoulaineEmpreintesInput,1),i,'rx');
%             
%             subplot(3,3,7);
%             hold on;
%             plot(D1.Convergence(1,:));
%             plot(D2.Convergence(1,:));
%             subplot(3,3,8);
%             hold on;
%             plot(D1.Convergence(2,:));
%             plot(D2.Convergence(2,:));
%             
%             i = i+2;
%     end
%     a = strcat(Path,'\Figures\');
%     cd(a)
%     savefig(strcat(num2str(i),'.fig'));
% end

for i = 1:S
    
    load(strcat(Path, Ressources(i,:)));
    D1 = Results;
    

%     D1.PoulaineEmpreintesInput(:,2) = -1* D1.PoulaineEmpreintesInput(:,2);
%     D1.PoulaineEmpreintesInput(:,5) = -1* D1.PoulaineEmpreintesInput(:,5);
%     D1.PoulaineEmpreintesInputScaled(:,2) = -1* D1.PoulaineEmpreintesInputScaled(:,2);
%     D1.PoulaineEmpreintesInputScaled(:,5) = -1* D1.PoulaineEmpreintesInputScaled(:,5);
%     D1.OriginalX(:,4) = -1* D1.OriginalX(:,4);
%     D1.OriginalX(:,5) = -1* D1.OriginalX(:,5);
    if Display
        DisplayCurves(D1.InitialPoulaine,i,'--');
        DisplayCurves(D1.PoulaineEmpreintesInput,i,'--');
        DisplayCurves(D1.PoulaineEmpreintesInputScaled,i,'--');
        DisplayCurves(D1.FinalPoulaine,i);
        DisplayCurves(D1.Mem.P,i);
        
        DisplayCurves(D1.Mem.TA,S+i);
        DisplayCurves(D1.TAFinal,S+i);
        DisplayCurves(D1.TAPostIK,S+i);
        
        %     DisplayCurves(D1.Convergence(1:3,:),3);
        
        
        sgtitle(Ressources(i,1:end-4),'Interpreter','none');
        %     IniTA = [];
        %     FinalTA = [];
        %
        
        
        IniTA = Pol2Curve(Results.InitialPolynom,size(Results.InitialPoulaine,1));
        FinalTA = Pol2Curve(Results.FinalPolynom,size(Results.InitialPoulaine,1));
        
        
        %     for j = 1:11
        %         IniTA = [IniTA, spline_to_curve_int(Results.InitialPolynom(:,1)==j,:),size(Results.InitialPoulaine,1)];
        %         FinalTA = [FinalTA , spline_to_curve_int(Results.FinalPolynom(Results.FinalPolynom(:,1)==j,:),size(Results.InitialPoulaine,1))];
        %     end
        
        f = figure(i);
        hold on;
        legend('Poulaine Initiale Scalée','Poulaine Empreintes','Poulaine Empreintes Scalée','Poulaine Finale','','','','','','','Position',[ 0.75, 0.1, 0.1, 0.1]);
        DisplayX(D1.OriginalX, size(D1.PoulaineEmpreintesInput,1),i,'ro');
        DisplayX(D1.EmpreintesScaled, size(D1.PoulaineEmpreintesInput,1),i,'rx');
        
        subplot(3,3,7);
        plot(D1.Convergence(1,:));
        subplot(3,3,8);
        plot(D1.Convergence(2,:));
        cd(strcat(Path,'\Figures\'))
        saveas(f,strcat(num2str(i),'.png'));
        
        %     g = figure(S+i);
        %     DisplayCurves(IniTA,S+i);
        %     DisplayCurves(FinalTA,S+i);
    end
end


Distances = struct;

A = [];
B = [];
C = [];
D = [];
E = [];
F = [];
G = [];
H = [];

for i = 1:S
    load(strcat(Path, Ressources(i,:)));
    
    Results.FinalPoulaine = Results.Mem.P;
    Results.TAFinal = Results.Mem.TA;
    
    if size(Results.FinalPoulaine,1) ~= size(Results.InitialPoulaine,1)
        Results.FinalPoulaine = Results.FinalPoulaine(1:end-1,:) ;
        Results.TAFinal = Results.TAFinal(1:end-1,:);
    end
    
    IniTA = spline_to_curve_int(Results.InitialPolynom,size(Results.InitialPoulaine,1));
%     FinalTA = spline_to_curve_int(Results.FinalPolynom,size(Results.InitialPoulaine,1));
    
    FinalTA = spline_to_curve_int(Results.Mem.Pol,size(Results.InitialPoulaine,1));
    
    A = [A, Results.Mem.Conv(1,end)];
    B = [B, Results.Mem.Conv(2,end) / Results.Convergence(2,1)];
    C = [C, sum(sum(sqrt((Results.FinalPoulaine - Results.InitialPoulaine).^2)))];
    D = [D, sum(sum(sqrt((Results.FinalPoulaine - Results.InitialSplinedPoulaine).^2)))];
    E = [E, Results.Convergence(1,1)];
%     F = [F, Results.Convergence(1,1)];
    % need to scale Input poulaine 2
    %         Distances(i).RToTargetScaledPoulaine = sqrt((Results.FinalPoulaine - Results.PoulaineEmpreintesInputScaled).^2);
    F = [F; sqrt(mean((-1* Results.FinalPoulaine + Results.InitialPoulaine).^2)) ];
    G = [G; sqrt(mean((-1* Results.FinalPoulaine + Results.InitialSplinedPoulaine).^2)) ];
    H = [H; sqrt(mean((-1* Results.TAFinal + Results.TAPostIK).^2))];    
end
% Distances.IToEmpreintes = F;
Distances.RToEmpreintes = A;
Distances.IniToEmpreintes = E;
Distances.DistCostVar = A./E;
Distances.MechCostChange = B;
Distances.RToOriginalPoulaine = C;
Distances.RToOriginalOptimPoulaine = D;
Distances.RMSP = F;
Distances.RMSPNS = G;
Distances.RMSTAIK =  H;
% Distances.DistRTAtoITA = sum(sum(sqrt((IniTA - FinalTA.^2))));

Stats.MeanDistToEmpreintes = mean(Distances.RToEmpreintes);
Stats.StdDistToEmpreintes = std(Distances.RToEmpreintes);
Stats.MeanMechCostChange = mean(Distances.MechCostChange);
Stats.StdMechCostChange = std(Distances.MechCostChange);
Stats.MeanDistToOriginalPoulaine = mean(Distances.RToOriginalPoulaine);
Stats.StdDistToOriginalPoulaine = std(Distances.RToOriginalPoulaine);
Stats.MeanDistToOriginalOptimPoulaine = mean(Distances.RToOriginalOptimPoulaine);
Stats.StdDistToOriginalOptimPoulaine = std(Distances.RToOriginalOptimPoulaine);
Stats.MeanDistCostVar = mean(Distances.DistCostVar);
Stats.StdDistCostVar = std(Distances.DistCostVar);
Stats.MeanDistIniToEmpreintes = mean(Distances.IniToEmpreintes);
Stats.StdDistIniToEmpreintes = std(Distances.IniToEmpreintes);

Stats.MeanRMSPoulaine = mean(Distances.RMSP);
Stats.StdRMSPoulaine = std(Distances.RMSP);
Stats.MeanRMSPoulaineNoSpline = mean(Distances.RMSPNS);
Stats.StdRMSPoulaineNoSpline = std(Distances.RMSPNS);
Stats.MeanRMSTA = mean(Distances.RMSTAIK);
Stats.StdRMSTA = std(Distances.RMSTAIK);

save(strcat(Path,'Stats.mat'),'Stats','Distances');

end

