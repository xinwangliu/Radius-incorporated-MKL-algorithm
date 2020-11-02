function [radius,radiusp,cost,Alpsupaux,w0aux,posaux] = costBasicRadiusclass(K0,StepSigma,DirSigma,Sigma,indsup,Alpsup,C,yapp,option)

global nbcall
nbcall=nbcall+1;

Sigma = Sigma+ StepSigma * DirSigma;
% kerneloption.matrix=sumKbeta(K,Sigma);
% kernel='numerical';

lambda=option.lambdareg;
verbose=option.verbosesvm;
alphainit=zeros(size(yapp));
alphainit(indsup)=yapp(indsup).*Alpsup;
[radius,radiusp,xsup,Alpsupaux,w0aux,posaux,alpha,cost] = basicRadiusclass([],yapp,C,lambda,K0,Sigma,verbose,alphainit);