function [cost,Alpsupaux,w0aux,posaux,nbsv] = costsvmoneagainstall(K,StepSigma,DirSigma,Sigma,C,yapp)

Sigma = Sigma+ StepSigma * DirSigma;
Kmatrix=sumKbeta(K,Sigma);

% warmstart.nbsv=nbsv;
% warmstart.alpsup=Alpsup;
% warmstart.pos=pos;
% verbose=option.verbosesvm;
% [xsup,Alpsupaux,w0aux,nbsv,posaux,cost]=svmmulticlassoneagainstall([],yapp,nbclass,C,lambdareg,kernel,kerneloption,verbose,warmstart);

[Alpsupaux,w0aux,nbsv,posaux,cost]=mySVMmulticlassoneagainstall(yapp,C,Kmatrix);
