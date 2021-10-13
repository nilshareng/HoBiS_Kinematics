function cost = KinforMin(Angles,Target,Param,Reperes)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
cost = norm(Target - fcineshort(Angles,Param,Reperes));
end

