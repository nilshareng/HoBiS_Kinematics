clear all;
close all;
clc;


XLPath = 'C:\Users\nhareng\Downloads\Poulaine_Fuite&Frimousse_treadmill.xls';

OutPath = 'C:\Users\nhareng\Desktop\CodeCommente\hobis\Ressources\PoulainesBabouin\';

[~, sheets] = xlsfinfo(XLPath);

Data = {};
Raw = struct;
for i = 1:size(sheets,2)
    Data{i} = xlsread(XLPath,sheets{i});
    
    [~,b] = size(Data{i});
    s = (b-4) / 3;
    Raw(i).X = [Data{i}(:,1:s),  Data{i}(:,end-2),] ;
    Raw(i).Y = [Data{i}(:,s+1:2*s),  Data{i}(:,end-1),] ;
    Raw(i).Z = [Data{i}(:,2*s+1:3*s),  Data{i}(:,end),] ;
    
    mkdir(strcat(OutPath,sheets{i}));
    
    for j = 1:s+1
        PN = [-1*Raw(i).Y(:,j), Raw(i).X(:,j), Raw(i).Z(:,j)]/10;
        PN = GaitSymetrisation(PN);
        dlmwrite(strcat(OutPath,sheets{i},'\',int2str(j),'.txt'),PN);
    end
    
end

