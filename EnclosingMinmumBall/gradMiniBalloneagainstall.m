function [grad] = gradMiniBalloneagainstall(K,pos,Alpsup,nbsv,radius,optimalBeta,Sigma)

nbclass=length(nbsv);
d = size(K,3);
nbsv=[0 nbsv];
aux=cumsum(nbsv);
S = zeros(1,d);
for i=1:nbclass
    waux=Alpsup(aux(i)+1:aux(i)+nbsv(i+1));
    indsup=pos(aux(i)+1:aux(i)+nbsv(i+1));
    S = S + gradMiniBallClass(K,indsup,waux,radius,optimalBeta,Sigma);
end
grad =S;
