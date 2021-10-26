function [] = DisplayGait(GaitMarkers,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    f=figure(1);
elseif nargin == 2
    f=figure(N);
end
hold on;
f.WindowState='maximized';

C = struct2cell(GaitMarkers);
field = fieldnames(GaitMarkers); 
for i = 1:size(C,3)
    tmp =  cell2struct(C(:,:,i),field,1); 
    
    if any(strcmp('RTarget',fieldnames(GaitMarkers)))
        plot3(tmp.RTarget(1),tmp.RTarget(2),tmp.RTarget(3),'rx');
    end
%     hold on;
    if any(strcmp('LTarget',fieldnames(tmp)))
        plot3(tmp.LTarget(1),tmp.LTarget(2),tmp.LTarget(3),'rx');
    end
    if any(strcmp('CRTarget',fieldnames(tmp)))
        plot3(tmp.CRTarget(1),tmp.CRTarget(2),tmp.CRTarget(3),'bo');
    end
    if any(strcmp('CLTarget',fieldnames(tmp)))
        plot3(tmp.CLTarget(1),tmp.CLTarget(2),tmp.CLTarget(3),'bo');
    end
    if any(strcmp('LPoul',fieldnames(tmp))) && any(strcmp('RPoul',fieldnames(tmp)))
        
        plot3(tmp.LPoul(:,1),tmp.LPoul(:,2),tmp.LPoul(:,3),'b');
%         f = figure('WindowState','maximized');
        xlim([-600 600]);
        ylim([-600 600]);
        zlim([-1000 200]);
        
        hold on;
        plot3(tmp.RPoul(:,1),tmp.RPoul(:,2),tmp.RPoul(:,3),'r');
    end

    P = DisplayModel(tmp,N);
    pause(1/20);
    hold off;
%     delete(P);
end

end

