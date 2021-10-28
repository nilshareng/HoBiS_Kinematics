function res = ErrorFun3(PolA,X,Sequence,Markers,Reperes)
% Renvoie la position de cheville associée aux TA aux instants de X selon les Données physiologiques Param 
tmpsum=[];
for i=1:size(X,1)
    %%% Eval des angles à tX
    tmpa=  zeros(1,11);
    for j = 1:11
        tmpa(j) = EvalSpline(PolA(PolA(:,1)==j,:),X(i,2));
    end
    
    
    tmp = fcine_numerique_H2(tmpa,Sequence,Markers,Reperes);
    if(X(i,1)<4)
        tmpsum = [tmpsum ; tmp(1:3)'];
    else
        tmpsum = [tmpsum ; tmp(4:6)'];
    end
end
res = tmpsum;
end


