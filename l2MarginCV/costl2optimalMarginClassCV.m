function [obj, Alpsupaux,w0aux,posaux] = costl2optimalMarginClassCV(K,StepSigma,DirSigma,Sigma,indsup,Alpsup,C,yapp,option)


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
[xsup,Alpsupaux,w0aux,posaux,alpha,obj] = l2optimalMarginClassCV([],yapp,C,lambdareg,kernel,kerneloption,verbosesvm,alphainit);
