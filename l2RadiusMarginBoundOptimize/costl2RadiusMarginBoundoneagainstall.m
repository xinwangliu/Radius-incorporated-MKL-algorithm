function [radius,optimalBeta,cost,Alpsupaux,w0aux,posaux,nbsv] = costl2RadiusMarginBoundoneagainstall(K,StepSigma,DirSigma,Sigma,Alpsup,yapp,pos,nbsv,nbclass,option)


Sigma = Sigma+ StepSigma * DirSigma;
lambdareg=option.lambdareg;
warmstart.nbsv=nbsv;
warmstart.alpsup=Alpsup;
warmstart.pos=pos;
verbose=option.verbosesvm;

[radius,optimalBeta,xsup,Alpsupaux,w0aux,nbsv,posaux,cost]= l2RadiusMarginBoundmulticlassoneagainstall([],yapp,nbclass,lambdareg,K, Sigma, verbose,warmstart);
