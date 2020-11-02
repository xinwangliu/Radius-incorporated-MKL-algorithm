function [cost, Alpsupaux,w0aux,posaux] = costl2MarginClass(K,StepSigma,DirSigma,Sigma,indsup,Alpsup,yapp,option)


global nbcall
nbcall=nbcall+1;

Sigma = Sigma+ StepSigma * DirSigma;


lambdareg=option.lambdareg;
verbosesvm=option.verbosesvm;
alphainit=zeros(size(yapp));
alphainit(indsup)=yapp(indsup).*Alpsup;

%compute margin
[xsup,Alpsupaux,w0aux,posaux,alpha,cost] = l2MarginClass([],yapp,lambdareg,K,Sigma,verbosesvm,alphainit);
