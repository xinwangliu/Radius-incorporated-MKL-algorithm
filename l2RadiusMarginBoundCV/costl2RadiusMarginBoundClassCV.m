function [obj,radius, optimalBeta, Alpsupaux,w0aux,posaux,margin] = costl2RadiusMarginBoundClassCV(K,StepSigma,DirSigma,Sigma,indsup,Alpsup,C,yapp,option)


global nbcall
nbcall=nbcall+1;

Sigma = Sigma+ StepSigma * DirSigma;
kerneloption.matrix=sumKbeta(K,Sigma);

kernel='numerical';
lambdareg=option.lambdareg;
verbosesvm=option.verbosesvm;
alphainit=zeros(size(yapp));
alphainit(indsup)=yapp(indsup).*Alpsup;

%compute margin
[radius, optimalBeta,xsup,Alpsupaux,w0aux,posaux,alpha,margin,obj] = l2MarginClassCV([],yapp,C,lambdareg,kernel,kerneloption,verbosesvm,alphainit);
