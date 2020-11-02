function [grad] = myGradsvmoneagainstall(K,pos,Alpsup,nbsv)

d=size(K,3);
grad = zeros(d,1);
nbclass=length(nbsv);
nbsv=[0 nbsv];
aux=cumsum(nbsv);
for k=1:d
    S = 0;
    for i=1:nbclass
        waux = Alpsup(aux(i)+1:aux(i)+nbsv(i+1));
        indsup=pos(aux(i)+1:aux(i)+nbsv(i+1));
        S=S + (- 0.5* waux'* K(indsup,indsup,k)*waux);
    end
    grad(k) = S;
end