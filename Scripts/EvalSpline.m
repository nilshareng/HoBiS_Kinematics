function res = EvalSpline(Pol, EvalPt)
% Evalue la valeur d'une spline � un point en t donn�  
% Pol array de type [ k PC1 PC2 Pol deg 3 ]
% EvalPt Pt d'�val

EvalPt = round(EvalPt,4);
% Pol(2:3,:) = round(Pol(2:3,:),4);
if(EvalPt<Pol(1,2))
    EvalPt= EvalPt+1-eps; %%%1-Periodicit�
end
% Si le polyn�me est bien conforme (non vide, un seul intervalle existe contenant EvalPt)
if (size(Pol(((Pol(:,2)<=EvalPt) & (Pol(:,3)>EvalPt)),4:7),1)~=1 || size(Pol(((Pol(:,2)<=EvalPt) & (Pol(:,3)>EvalPt)),4:7),2)~=4)
    res= -10;
    return
end
% Evaluation du polynome
res = polyval(Pol(((Pol(:,2)<=EvalPt) & (Pol(:,3)>EvalPt)),4:7),EvalPt);

end

