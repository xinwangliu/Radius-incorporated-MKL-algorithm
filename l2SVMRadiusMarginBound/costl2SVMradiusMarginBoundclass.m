function [cost,radius, optimalBeta, Alpsupaux,w0aux,posaux, margin] = costl2SVMradiusMarginBoundclass(K,StepSigma,DirSigma,Sigma,indsup,Alpsup,C,yapp,option)


global nbcall
nbcall=nbcall+1;

Sigma = Sigma+ StepSigma * DirSigma;
kerneloption.matrix=sumKbeta(K,Sigma);
lambda=1e-10;
verbose=0;

%compute Radius;
[radius, optimalBeta]= l2SVMradius(kerneloption.matrix,C,lambda,verbose);

kernel='numerical';
span=1;
lambdareg=option.lambdareg;
verbosesvm=option.verbosesvm;
alphainit=zeros(size(yapp));
alphainit(indsup)=yapp(indsup).*Alpsup;

%compute margin
[xsup,Alpsupaux,w0aux,posaux,alpha,margin] = marginClass([],yapp,C,lambdareg,kernel,kerneloption,verbosesvm,alphainit);

% compute the cost
cost = radius * margin;