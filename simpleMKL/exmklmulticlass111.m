%
% Example MKL MultiClass SVM Classifiction
%

close all
clear all
clc
%------------------------------------------------------
% Creating data
%------------------------------------------------------
n=20;
sigma=1.2;
nbclass=3;

x1= sigma*randn(n,2)+ ones(n,1)*[-1.5 -1.5];
x2= sigma*randn(n,2)+ ones(n,1)*[0 2];
x3= sigma*randn(n,2)+ ones(n,1)*[2 -1.5];
xapp=[x1;x2;x3];
yapp=[1*ones(1,n) 2*ones(1,n) 3*ones(1,n)]';

[n1, n2]=size(xapp);
[xtesta1,xtesta2]=meshgrid([-4:0.1:4],[-4:0.1:4]);
[na,nb]=size(xtesta1);
xtest1=reshape(xtesta1,1,na*nb);
xtest2=reshape(xtesta2,1,na*nb);
xtest=[xtest1;xtest2]';

addpath('/home/users/xliu/XINWAN_WORK/MKLToolBox/mkLearningWithRadius/');
addpath('/home/users/xliu/XINWAN_WORK/MKLToolBox/SVM-KM/');
addpath('/home/users/xliu/XINWAN_WORK/MKLToolBox/simpleMKL/');
%----------------------------------------------------------
%   Learning and Learning Parameters
%   Parameters are similar to those used for mklsvm
%-----------------------------------------------------------

C = 100;
lambda = 1e-7;
verbose = 1;
options.algo='oneagainstall';
options.seuildiffsigma=1e-4;
options.seuildiffconstraint=0.1;
options.seuildualitygap=1e-2;
options.goldensearch_deltmax=1e-1;
options.numericalprecision=1e-8;
options.stopvariation=1;
options.stopKKT=0;
options.stopdualitygap=0;
options.firstbasevariable='first';
options.nbitermax=500;
options.seuil=0.;
options.seuilitermax=10;
options.lambdareg = 1e-6;
options.miniter=0;
options.verbosesvm=0;
options.efficientkernel=0;
%------------------------------------------------------------

kernelt={'gaussian' 'gaussian' 'poly' 'poly' };
kerneloptionvect={[0.5 1 2 5 7 10 12 15 17 20] [0.5 1 2 5 7 10 12 15 17 20] [1 2 3] [1 2 3]};
variablevec={'all' 'single' 'all' 'single'};



[nbdata,dim]=size(xapp);
[kernel,kerneloptionvec,variableveccell]=CreateKernelListWithVariable(variablevec,dim,kernelt,kerneloptionvect);
[Weight,InfoKernel]=UnitTraceNormalization(xapp,kernel,kerneloptionvec,variableveccell);
K=mklkernel(xapp,InfoKernel,Weight,options);
nbkernel = size(K,3);

%---------------------Learning & Testing ----------------


[beta,w,w0,pos,nbsv,SigmaH,obj] = mklmulticlass(K,yapp,C,nbclass,options,verbose);
xsup=xapp(pos,:);
Kt=mklkernel(xtest,InfoKernel,Weight,options,xsup,beta);
kernel='numerical';
kerneloption.matrix=Kt;

[ypred,maxi] = svmmultival([],[],w,w0,nbsv,kernel,kerneloption);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K0 = zeros(nbdata,nbdata,nbkernel+1);
K0(:,:,1:nbkernel) = K;
K0(:,:,nbkernel+1) = eye(nbdata)/nbdata;

[Sigma1,Alpsup1,w01,pos1,nbsv1,SigmaH1,obj1] = l2trStmklmulticlass(K0,yapp,nbclass,options,verbose);
xsup1=xapp(pos1,:);
Kt1=mklkernel(xtest,InfoKernel,Weight,options,xsup1,Sigma1(1:end-1));
kerneloption1.matrix=Kt1;

[ypred1,maxi1] = svmmultival([],[],Alpsup1,w01,nbsv1,kernel,kerneloption1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[radius,Sigma2,Alpsup2,w02,pos2,nbsv2,SigmaH2,obj2] = miniBallmklmulticlass(K,yapp,C,nbclass,options,verbose);
xsup2=xapp(pos2,:);
Kt2=mklkernel(xtest,InfoKernel,Weight,options,xsup2,Sigma2);
kerneloption2.matrix=Kt2;

[ypred2,maxi2] = svmmultival([],[],Alpsup2/radius,w02,nbsv2,kernel,kerneloption2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Sigma4,Alpsup4,w04,pos4,nbsv4,SigmaH4,obj4] = ecml2009mklmulticlass(K,yapp,C,nbclass,options,verbose);
% xsup4=xapp(pos4,:);
% Kt4=mklkernel(xtest,InfoKernel,Weight,options,xsup4,Sigma4);
% kerneloption4.matrix=Kt4;
% 
% [ypred4,maxi4] = svmmultival([],[],Alpsup4,w04,nbsv4,kernel,kerneloption4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[radius3,Sigma3,Alpsup3,w03,pos3,nbsv3] = l2RadiusMarginBoundmklmulticlass(K0,yapp,nbclass,options,verbose);
xsup3=xapp(pos3,:);
Kt3=mklkernel(xtest,InfoKernel,Weight,options,xsup3,Sigma3(1:end-1));
kerneloption3.matrix=Kt3;

[ypred3,maxi3] = svmmultival([],[],Alpsup3/radius3,w03,nbsv3,kernel,kerneloption3);


% %------------------------------------------------------------------
% %           Plotting the decision function
% %-------------------------------------------------------------------
% ypredmat=reshape(ypred,na,nb);
% contour(xtesta1,xtesta2,ypredmat,[1 2 3]);hold on
% style=['x+*'];
% color=['bgr'];
% hold on
% for i=0:nbclass-1
%     h=plot(xapp(i*n+1:(i+1)*n,1),xapp(i*n+1:(i+1)*n,2),[style(i+1) color(i+1)]);
%     set(h,'LineWidth',2);
%     hold on
% end;

% if ~isempty(xsup)
%     h=plot(xsup(:,1),xsup(:,2),'ok');
%     set(h,'LineWidth',2);
%     axis( [ -4 4 -4 4]);
%     legend('classe 1','classe 2','classe 3', 'Support Vector');
%     hold off
% end




