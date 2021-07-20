%%
Manipulability = struct;


% for i = 1:(floor(size(GT,1)))
%  
%     
%     
%     [Jk, V, DJerk] = calc_jacobien_PC_3D(NPCA, NPolA,derflag, 0, X, dp,M, Cost, JerkRef, Fem1g, Fem1d, Fem6g, Fem6d, Tal1g, Tal1d, R_monde_local,R_Pelvis_monde_local, R_LFem_ref_local, R_LTib_ref_local, R_RFem_ref_local, R_RTib_ref_local);
%     
%     Manipulability(i).Jx=Jk(1:6,:);
%     Manipulability(i).Jy=Jk(7:12,:);
%     Manipulability(i).Jz=Jk(13:end,:);
% 
% end

% 
% for i =1:size(Manipulability,2)
%     Manipulability(i).Kinx = kmanip(Manipulability(i).Jx);
% end

Data = struct;

Data.Param = Param;
Data.Reperes= Reperes;
% Data.Repere.R_monde_local = R_monde_local ;
% Data.Repere.R_Pelvis_monde_local = R_Pelvis_monde_local;
% Data.Repere.R_LFem_ref_local = R_LFem_ref_local;
% Data.Repere.R_LTib_ref_local = R_LTib_ref_local;
% Data.Repere.R_RFem_ref_local = R_RFem_ref_local;
% Data.Repere.R_RTib_ref_local = R_RTib_ref_local;

Data.X = X;
Data.PCA = NPCA;
Data.PolA = NPolA;
Data.Poulaine = NewPoul;
Data.NewCurve = NewCurve;
Data.Manip = Manipulability;

PolP = [PolP(:,1) , PolP(:,3:end) ];

Data.PolP = PolP;

for i  = 1:size(GT,1)
    Data.EvalTime = (i-1)/(size(GT,1)-1);
    Manipulability(i).J = JacMan(Data);
    Manipulability(i).Kin = kmanip(Manipulability(i).J);
end

Data.Manip = Manipulability;


tmp=[];
SpeedCurve=[];
NPAxis=[];



for i =1:size(Manipulability,2)
    NPAxis(:,i) = Manipulability(i).Kin.PAxis'/norm( Manipulability(i).Kin.PAxis);
    dPolP = PolP;
    dPolP(:,4) = zeros(size(dPolP(:,4)));
    dPolP(:,5) = 3*PolP(:,4);
    dPolP(:,6) = 2*PolP(:,5);
    dPolP(:,7) = PolP(:,6);
    
    tmp=[];
    for j=1:6
        tmp = [tmp , EvalSpline(dPolP(dPolP(:,1)==j,:),(i-1)/(size(Manipulability,2)-1))];
    end
    SpeedCurve = [SpeedCurve ; tmp];
end

Data.dPolP = dPolP;
Data.SpeedCurve = SpeedCurve;

save(strcat('C:\Users\nhareng\Desktop\CodeCommente\ModelePied\ManipRes2\Manip',int2str(jj),'.mat'),'Data');

% %%
% 
% clear all; close all; clc;
% 
% load('C:\Users\nhareng\Desktop\CodeCommente\ModelePied\Manip2.mat')
% 

% 
% % for i =1:Period+1
% %     tmp=zeros(1,11);
% %     for j=1:11
% %         tmp(j) = EvalSpline(Pol(Pol(:,1)==j,:),((i-1))/Period);
% %     end
% %     NewCurve = [NewCurve ; tmp];
% % end
% 
% 
% % for i = 1:size(Manipulability,2)
% %     Anim3DParam(NewCurve,n,v,Param)
% % end
% % 
% % figure;

