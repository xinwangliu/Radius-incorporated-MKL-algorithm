function [Kx] = myCalculateKx(w,pos,nbsv,num)

Kx = zeros(num);
nbclass = length(nbsv);
nbsv =[0 nbsv];
aux = cumsum(nbsv);
for ic = 1: nbclass
    alphaic = zeros(num,1);
    waux = w(aux(ic)+1:aux(ic)+nbsv(ic+1));
    posi = pos(aux(ic)+1:aux(ic)+nbsv(ic+1));
    alphaic(posi) = waux;
    Kx = Kx + alphaic*alphaic';
end