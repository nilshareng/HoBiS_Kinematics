function [Gait, GaitMarkers, GaitReperes] = Angles2Gait(Angles,Sequence,Markers,Reperes)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
GaitReperes = [];
GaitMarkers = [];
Gait = zeros(size(Angles,1),6);

for i = 1:size(Angles,1)
   [Gait(i,:), TmpMarkers , TmpReperes] = fcinematique(Angles(i,:),Sequence,Markers,Reperes);
    GaitMarkers = [GaitMarkers, TmpMarkers];
    GaitReperes = [GaitReperes, TmpReperes];
end
end

