function [cost,radius, optimalBeta, Alpsupaux,w0aux,posaux, margin] = costRadiusMarginBoundClass(K,StepSigma,DirSigma,Sigma,indsup,Alpsup,yapp,option)


global nbcall
nbcall=nbcall+1;

Sigma = Sigma+ StepSigma * DirSigma;


lambdareg=option.lambdareg;
verbosesvm=option.verbosesvm;
alphainit=zeros(size(yapp));
alphainit(indsup)=yapp(indsup).*Alpsup;

%compute margin and radius
[radius, optimalBeta,xsup,Alpsupaux,w0aux,posaux,alpha,margin,cost] = radiusMarginBoundClass([],yapp,lambdareg,K,Sigma,verbosesvm,alphainit);

