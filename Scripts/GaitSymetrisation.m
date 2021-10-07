function SymCurves = GaitSymetrisation(Curves)
% Symmetricalisation of gait curves, for ankles traj and angular traj
% Ankle traj : 2*3 DoF 
% Angular Traj : 11 DoF, 3 + 2*4 
SymCurves = [];
switch size(Curves,2)
    
    case {3, 6} % Ankles, phase opposition for X, Y and Z ; Amp opposition for X 
        mid = fix(size(Curves,1)/2);
        SymCurves = [Curves(:,1:3), Curves(:,1:3)];
        SymCurves(:,4) = -1*[SymCurves(mid+1:end,1) ; SymCurves(1:mid,1)];
        SymCurves(:,5) = [SymCurves(mid+1:end,2) ; SymCurves(1:mid,2)];
        SymCurves(:,6) = [SymCurves(mid+1:end,3) ; SymCurves(1:mid,3)];
                        
    case {7, 11}
        mid = fix(size(Curves,1)/2);
        SymCurves = [Curves(:,1:7) , Curves(:,4:7)];
        SymCurves(:,8) = -1*[SymCurves(mid+1:end,8) ; SymCurves(1:mid,8)];
        SymCurves(:,9) = -1*[SymCurves(mid+1:end,9) ; SymCurves(1:mid,9)];
        SymCurves(:,10) = [SymCurves(mid+1:end,10) ; SymCurves(1:mid,10)];
        SymCurves(:,11) = [SymCurves(mid+1:end,11) ; SymCurves(1:mid,11)];
end

