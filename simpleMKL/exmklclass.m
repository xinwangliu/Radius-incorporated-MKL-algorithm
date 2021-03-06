% Example of how to use the mklsvm for  classification
clear 
clc

addpath(genpath('C:\XINWANGWORK\MachineLearningToolbox\mkLearningWithRadius\KernelCode\simpleMKL\'));
addpath(genpath('C:\XINWANGWORK\MachineLearningToolbox\MaifoldRegMKL\'));
addpath(genpath('C:\XINWANGWORK\MachineLearningToolbox\hessmkl\'));

nbiter=1;
ratio=0.05;
data='ionosphere';
verbose=1;
options.algo='svmclass'; % Choice of algorithm in mklsvm can be either
                         % 'svmclass' or 'svmreg'
%------------------------------------------------------
% choosing the stopping criterion
%------------------------------------------------------
options.stopvariation=1; % use variation of weights for stopping criterion 
options.stopKKT=0;       % set to 1 if you use KKTcondition for stopping criterion    
options.stopdualitygap=0; % set to 1 for using duality gap for stopping criterion

%------------------------------------------------------
% choosing the stopping criterion value
%------------------------------------------------------
options.seuildiffsigma=1e-4;        % stopping criterion for weight variation 
options.seuildiffconstraint=0.1;    % stopping criterion for KKT
options.seuildualitygap=0.01;       % stopping criterion for duality gap

%------------------------------------------------------
% Setting some numerical parameters 
%------------------------------------------------------
options.goldensearch_deltmax=1e-1; % initial precision of golden section search
options.numericalprecision=1e-16;   % numerical precision weights below this value
                                   % are set to zero 
options.lambdareg = 1e-8;          % ridge added to kernel matrix 

%------------------------------------------------------
% some algorithms paramaters
%------------------------------------------------------
options.firstbasevariable='first'; % tie breaking method for choosing the base 
                                   % variable in the reduced gradient method 
options.nbitermax=500;             % maximal number of iteration  
options.seuil=0;                   % forcing to zero weights lower than this 
options.seuilitermax=10;           % value, for iterations lower than this one 

options.miniter=0;                 % minimal number of iterations 
options.verbosesvm=0;              % verbosity of inner svm algorithm 
options.efficientkernel=0;         % use efficient storage of kernels 
options.threshold = 1e-4;

%------------------------------------------------------------------------
%                   Building the kernels parameters
%------------------------------------------------------------------------
kernelt={'gaussian' 'poly' }; %%{'gaussian' 'gaussian' 'poly' 'poly' };
kerneloptionvect={[0.5 1 2 5 7 10 12 15 17 20] [1:5]};
variablevec={'all' 'all'};  %%{'all' 'single' 'all' 'single'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classcode=[1 -1];
load([data ]);
[nbdata,dim]=size(x);

nbtrain=floor(nbdata*ratio);
rand('state',0);

rho = 100;
C = 10;
for num = 1: 1
    num
    [xapp,yapp,xtest,ytest,indice]=CreateDataAppTest(x, y, nbtrain,classcode);
    [xapp,xtest]=normalizemeanstd(xapp,xtest);
    [kernel,kerneloptionvec,variableveccell]=CreateKernelListWithVariable(variablevec,dim,kernelt,kerneloptionvect);
    %%%KTrn
    [Weight,InfoKernel]=UnitTraceNormalization(xapp,kernel,kerneloptionvec,variableveccell);
    K = mklkernel(xapp,InfoKernel,Weight,options);
    KW = mklkernel([xapp;xtest],InfoKernel,Weight,options);
    nbkernel = size(K,3);
  
    %%%%Laplacian Matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% SimpleMKL
    tic
    [Sigma1,Alpsup1,w01,pos1,obj1] = mklsvm(K,yapp,C,options,verbose);
    timelasso1(num)=toc;

    Kt1=mklkernel(xtest,InfoKernel,Weight,options,xapp(pos1,:),Sigma1);
    ypred1=Kt1*Alpsup1 + w01;
    bc1(num)=mean(sign(ypred1)==ytest);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options.qnorm = 2;
    tic
    [Sigma2,Alpsup2,w02,pos2,obj2] = mklGLasso_alt0_class(K,yapp,C,options,verbose);
    timelasso2(num)=toc;
    
    Kt2=mklkernel(xtest,InfoKernel,Weight,options,xapp(pos2,:),Sigma2);
    ypred2=Kt2*Alpsup2 + w02;
    bc2(num)=mean(sign(ypred2)==ytest);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    options.qnorm = 3;
    tic
    [Sigma3,Alpsup3,w03,pos3,obj3] = mklGLasso_alt0_class(K,yapp,C,options,verbose);
    timelasso3(num)=toc;
    
    Kt3=mklkernel(xtest,InfoKernel,Weight,options,xapp(pos3,:),Sigma3);
    ypred3=Kt3*Alpsup3 + w03;
    bc3(num)=mean(sign(ypred3)==ytest);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    [L] = myLaplacian(KW);
    Sigma4 = kernel_combination(K,L,yapp,C,rho);
    [Alpsup4,w04,pos4]=hessianMKLclass(K,Sigma4,yapp,C);
    
    Kt4=mklkernel(xtest,InfoKernel,Weight,options,xapp(pos4,:),Sigma4);
    ypred4=Kt4*Alpsup4 + w04;
    bc4(num)=mean(sign(ypred4)==ytest);
 
end%



