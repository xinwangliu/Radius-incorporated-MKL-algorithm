function [ypred,maxi,ypredMat]= l2trStmultival(w,b,pos,nbsv,Ktsttrn,Sigma)

% USAGE ypred=svmmultival(x,xsup,w,b,nbsv,kernel,kerneloption)
%
% Process the class of a new point x of a one-against-all
% or a all data at once MultiClass SVM
%
Kmatrix = sumKbeta(Ktsttrn,Sigma);

nbclass=length(nbsv);
nbsv=[0 nbsv];
aux=cumsum(nbsv);
for i=1:nbclass
    %  Kernel matrix is given as a parameter
    waux=w(aux(i)+1:aux(i)+nbsv(i+1));
    baux=b(i);
    posi = pos(aux(i)+1:aux(i)+nbsv(i+1));
    %% [y] = mysvmval(w,b,pos,Ktsttrn)
    ypred(:,i)= mysvmval(waux,baux,posi,Kmatrix);
end;
ypredMat=ypred;
[maxi,ypred]=max(ypred');
maxi=maxi';
ypred=ypred';