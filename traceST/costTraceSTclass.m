function [trST,trSTp,cost,Alpsupaux,w0aux,posaux] = costTraceSTclass(K0,StepSigma,DirSigma,Sigma,indsup,Alpsup,C,yapp,option)

global nbcall
nbcall=nbcall+1;

Sigma = Sigma+ StepSigma * DirSigma;
% kerneloption.matrix=sumKbeta(K,Sigma);
% kernel='numerical';

lambda=option.lambdareg;
verbose=option.verbosesvm;
alphainit=zeros(size(yapp));
alphainit(indsup)=yapp(indsup).*Alpsup;
[trST,trSTp,xsup,Alpsupaux,w0aux,posaux,alpha,cost] = traceSTclass([],yapp,C,lambda,K0,Sigma,verbose,alphainit);