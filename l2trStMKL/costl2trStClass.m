function [cost,Alpsupaux,w0aux,posaux] = costl2trStClass(K0,StepSigma,DirSigma,Sigma,indsup,Alpsup,yapp,option)

global nbcall
nbcall=nbcall+1;

Sigma = Sigma+ StepSigma * DirSigma;

lambda=option.lambdareg;
verbose=option.verbosesvm;
alphainit=zeros(size(yapp));
alphainit(indsup)=yapp(indsup).*Alpsup;
[xsup,Alpsupaux,w0aux,posaux,alpha,cost] = l2trStClass([],yapp,lambda,K0,Sigma,verbose,alphainit);