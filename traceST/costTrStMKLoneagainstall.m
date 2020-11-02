function [trST,trSTp,cost,Alpsupaux,w0aux,posaux,nbsv] = costTrStMKLoneagainstall(K,C,StepSigma,DirSigma,Sigma,Alpsup,yapp,pos,nbsv,nbclass,option)


Sigma = Sigma+ StepSigma * DirSigma;
lambdareg=option.lambdareg;
warmstart.nbsv=nbsv;
warmstart.alpsup=Alpsup;
warmstart.pos=pos;
verbose=option.verbosesvm;

[trST,trSTp,xsup,Alpsupaux,w0aux,nbsv,posaux,cost]= trStmulticlassoneagainstall([],yapp,C,nbclass,lambdareg,K, Sigma, verbose,warmstart);
