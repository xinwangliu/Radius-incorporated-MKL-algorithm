function [cost,Alpsupaux,w0aux,posaux,nbsv] = costl2trStMKLoneagainstall(K,StepSigma,DirSigma,Sigma,Alpsup,yapp,pos,nbsv,nbclass,option)


Sigma = Sigma+ StepSigma * DirSigma;

lambdareg=option.lambdareg;

warmstart.nbsv=nbsv;
warmstart.alpsup=Alpsup;
warmstart.pos=pos;
verbose=option.verbosesvm;

[xsup,Alpsupaux,w0aux,nbsv,posaux,cost]=l2trStmulticlassoneagainstall([],yapp,nbclass,lambdareg,K, Sigma, verbose,warmstart);

