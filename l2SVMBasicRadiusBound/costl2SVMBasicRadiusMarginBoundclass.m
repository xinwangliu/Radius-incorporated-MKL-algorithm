function [cost,radius, radiusp, Alpsupaux,w0aux,posaux] = costl2SVMBasicRadiusMarginBoundclass(K,StepSigma,DirSigma,Sigma,indsup,Alpsup,C,yapp,option)


global nbcall
nbcall=nbcall+1;

Sigma = Sigma+ StepSigma * DirSigma;

lambdareg=option.lambdareg;
verbosesvm=option.verbosesvm;
alphainit=zeros(size(yapp));
alphainit(indsup)=yapp(indsup).*Alpsup;

%compute margin
[radius,radiusp, xsup,Alpsupaux,w0aux,posaux,alpha,cost] = l2SVMBasicRadiusClass([],yapp,C,lambdareg,K,Sigma,verbosesvm,alphainit);