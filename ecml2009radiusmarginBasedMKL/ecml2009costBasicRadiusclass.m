function [cost,Alpsupaux,w0aux,posaux] = ecml2009costBasicRadiusclass(K0,StepSigma,DirSigma,Sigma,C,radiusp,yapp,option)

global nbcall
nbcall=nbcall+1;

Sigma = Sigma+ StepSigma * DirSigma;

[Alpsupaux,w0aux,posaux,cost] = ecml2009basicRadiusclass(K0,Sigma,yapp,C,radiusp,option);
