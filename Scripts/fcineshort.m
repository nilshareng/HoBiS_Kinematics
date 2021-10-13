function [res,Save,markers] = fcineshort(Angles, Param, Reperes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if size(Param,2)~=3
    Param = Param';
end

if (exist markers)
    [res,Save,markers] = fcine_numerique_H3_markers(Angles, markers,Param(1,:),Param(2,:),Param(3,:),Param(4,:),Param(5,:),Param(6,:), Reperes.Monde, Reperes.Pelvis, Reperes.LFemurLocal, Reperes.LTibiaLocal, Reperes.RFemurLocal, Reperes.RTibiaLocal);
    
else
    
    %     [res,Save] = fcine_numerique_H3(Angles,Param(1,:),Param(2,:),Param(3,:),Param(4,:),Param(5,:),Param(6,:), Reperes.MondeLocal(1:3,1:3), Reperes.PelvisLocal(1:3,1:3), Reperes.LFemurLocal(1:3,1:3), Reperes.LTibiaLocal(1:3,1:3), Reperes.RFemurLocal(1:3,1:3), Reperes.RTibiaLocal(1:3,1:3));
    [res,Save] = fcine_numerique_H3(Angles,Param(1,:),Param(2,:),Param(3,:),Param(4,:),Param(5,:),Param(6,:), Reperes.Monde, Reperes.Pelvis, Reperes.LFemurLocal, Reperes.LTibiaLocal, Reperes.RFemurLocal, Reperes.RTibiaLocal);
end

end
