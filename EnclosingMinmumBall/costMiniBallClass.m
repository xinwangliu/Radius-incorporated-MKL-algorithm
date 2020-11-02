function [radius,optimalBeta,cost,Alpsupaux,w0aux,posaux] = costMiniBallClass(K,StepSigma,DirSigma,Sigma,indsup,Alpsup,C,yapp,option)

global nbcall
nbcall=nbcall+1;

Sigma = Sigma+ StepSigma * DirSigma;

lambda=option.lambdareg;
verbose=option.verbosesvm;
alphainit=zeros(size(yapp));
alphainit(indsup)=yapp(indsup).*Alpsup;
[radius,optimalBeta,xsup,Alpsupaux,w0aux,posaux,alpha,cost] = miniBallClass([],yapp,C,lambda,K,Sigma,verbose,alphainit);
