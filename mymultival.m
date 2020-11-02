function [ypred,label,maxi]= mymultival(w,b,pos,nbsv,Ktsttrn,Sigma)

% USAGE ypred=svmmultival(x,xsup,w,b,nbsv,kernel,kerneloption)
%
% Process the class of a new point x of a one-against-all
% or a all data at once MultiClass SVM
%
Kavertsttrn = sumKbeta(Ktsttrn,Sigma);

nbclass=length(nbsv);
nbsv=[0 nbsv];
aux=cumsum(nbsv);
ntest = size(Ktsttrn,1);
ypred = zeros(ntest,nbclass);
for i=1:nbclass
    %  Kernel matrix is given as a parameter
    waux=w(aux(i)+1:aux(i)+nbsv(i+1));
    baux=b(i);
    posi = pos(aux(i)+1:aux(i)+nbsv(i+1));
    %%% [y] = mysvmval(w,b,pos,Ktsttrn)
    ypred(:,i)= Kavertsttrn(:,posi)*waux+baux;
end
[maxi,label]=max(ypred,[],2);