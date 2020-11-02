function [radius,radiusp,cost,Alpsupaux,w0aux,posaux,nbsv] = costecml2009MKLoneagainstall(K,C,StepSigma,DirSigma,Sigma,Alpsup,yapp,pos,nbsv,nbclass,option)


Sigma = Sigma+ StepSigma * DirSigma;
warmstart.nbsv=nbsv;
warmstart.alpsup=Alpsup;
warmstart.pos=pos;

[radius,radiusp,xsup,Alpsupaux,w0aux,nbsv,posaux,cost]=  ecml2009multiclassoneagainstall([],yapp,C,nbclass,K, Sigma, option,warmstart);
