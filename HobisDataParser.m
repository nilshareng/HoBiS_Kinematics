function [markers, Param]= HobisDataParser(Data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
markers = struct;
dd = Data;

Param = zeros(6,3);
% Ancien parser de 'JeuEssaiDM'
Pelv1=dd(1,:);	
Pelv2g=dd(2,:);
Pelv4g=dd(3,:);
Pelv5g=dd(4,:);

Fem1g=dd(5,:);
Fem9g=dd(6,:);
Fem10g=dd(7,:);
Fem6g=dd(8,:);

Tib5g=dd(9,:);
Tib6g=dd(10,:);
Tib1g=dd(11,:);

Tal2g=dd(12,:);
Tal3g=dd(13,:);
Tal1g=dd(14,:);

Pelv2d=dd(15,:);
Pelv4d=dd(16,:);
Pelv5d=dd(17,:);

Fem1d=dd(18,:);
Fem9d=dd(19,:);
Fem10d=dd(20,:);
Fem6d=dd(21,:);

Tib5d=dd(22,:);
Tib6d=dd(23,:);
Tib1d=dd(24,:);

Tal2d=dd(25,:);
Tal3d=dd(26,:);
Tal1d=dd(27,:);

Param = [Fem1g',Fem6g',Tal1g', Fem1d', Fem6d', Tal1d']/1000;


markers.RFWT = Pelv4d;
markers.LFWT = Pelv4g;
markers.RBWT = Pelv5d;
markers.LBWT = Pelv5g;
markers.RKNE = Fem10d;
markers.LKNE = Fem10g;
markers.RKNI = Tib5d;
markers.LKNI = Tib5g;
markers.RANE = Tal2d;
markers.LANE = Tal2g;
markers.RANI = Tal3d;
markers.LANI = Tal3g;

end

