function [radius,optimalBeta,xsup,w,b,nbsv,pos,obj] = miniBallmulticlassoneagainstall(x,y,C,nbclass,lambda,K,Sigma,verbose,warmstart)

%  [xsup,w,d,pos,alpha,obj]=l2trStClass(x,y,lambda, K,Sigma,verbose,alphainit)

% SVM Multi Classes Classification One against Others algorithm
%
% y is the target vector which contains integer from 1 to nbclass.
% 
% This subroutine use the svmclass function
% 
% the differences lies in the output nbsv which is a vector
% containing the number of support vector for each machine
% learning.
% For xsup, w, b, the output of each ML are concatenated
% in line.
% 
%
% See svmclass, svmmultival
%
%
% 29/07/2000 Alain Rakotomamonjy


xsup=[];  % 3D matrices can not be used as numebr of SV changes
w=[];
b=[];
pos=[];
nbsv=zeros(1,nbclass);
obj=0;

for i=1:nbclass

    yone=(y==i)+(y~=i)*-1;
    if exist('warmstart') & isfield(warmstart,'nbsv');
        nbsvinit=cumsum([0 warmstart.nbsv]);
        alphainit=zeros(size(yone));
        alphainit(warmstart.pos(nbsvinit(i)+1:nbsvinit(i+1)))= abs(warmstart.alpsup(nbsvinit(i)+1:nbsvinit(i+1)));
    else
        alphainit=[];
    end;

    [radius, optimalBeta, xsupaux,waux,baux,posaux,alphaaux,objaux]= miniBallClass(x,yone,C,lambda,K,Sigma,verbose,alphainit);

    nbsv(i)=length(posaux);
    xsup=[xsup;xsupaux];
    w=[w;waux];
    b=[b;baux];
    pos=[pos;posaux];
    obj=obj+objaux;
end;
