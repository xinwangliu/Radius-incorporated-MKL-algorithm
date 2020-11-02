% Example of how to use the mklsvm for  classification
% clear all
% close all
% addpath('../toollp');
clear;
clc;
warning off;

nbiter=1;
ratio=0.1;
data='ionosphere';
C = [100];
verbose=1;
lambda0 = 2^5;
numClusterPos = 2;
numClusterNeg= 2;

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
options.seuildiffsigma=1e-3;        % stopping criterion for weight variation 
options.seuildiffconstraint=0.1;    % stopping criterion for KKT
options.seuildualitygap=0.01;       % stopping criterion for duality gap

%------------------------------------------------------
% Setting some numerical parameters 
%------------------------------------------------------
options.goldensearch_deltmax=1e-1; % initial precision of golden section search
options.numericalprecision=1e-8;   % numerical precision weights below this value
                                   % are set to zero 
options.lambdareg = 1e-8;          % ridge added to kernel matrix 

%------------------------------------------------------
% some algorithms paramaters
%------------------------------------------------------
options.firstbasevariable='first'; % tie breaking method for choosing the base 
                                   % variable in the reduced gradient method 
options.nbitermax=100;             % maximal number of iteration  
options.seuil=0;                   % forcing to zero weights lower than this 
options.seuilitermax=10;           % value, for iterations lower than this one 

options.miniter=0;                 % minimal number of iterations 
options.verbosesvm=0;              % verbosity of inner svm algorithm 
options.efficientkernel=0;         % use efficient storage of kernels 
%------------------------------------------------------------------------
%                   Building the kernels parameters
%------------------------------------------------------------------------
kernelt={'gaussian'};
kerneloptionvect={[2 4 5 6 7 8]};  %[0.5 1 2 5 7 10 12 15 17 20]
variablevec={'all'};

classcode=[1 -1];
load([data ]);
[nbdata,dim]=size(x);

nbtrain=floor(nbdata*ratio);
rand('state',0);

for i=1: nbiter
    i
    [xapp,yapp,xtest,ytest,indice]=CreateDataAppTest(x, y, nbtrain,classcode);
    [xapp,xtest]=normalizemeanstd(xapp,xtest);
    [kernel,kerneloptionvec,variableveccell]=CreateKernelListWithVariable(variablevec,dim,kernelt,kerneloptionvect);
    [Weight,InfoKernel]=UnitTraceNormalization(xapp,kernel,kerneloptionvec,variableveccell);
    K=mklkernel(xapp,InfoKernel,Weight,options);

    %------------------------------------------------------------------
    % 
    %  K is a 3-D matrix, where K(:,:,i)= i-th Gram matrix 
    %
    %------------------------------------------------------------------
    % or K can be a structure with uses a more efficient way of storing
    % the gram matrices
    %
    % K = build_efficientK(K);
    
    tic
    [Sigma,w,b,posw,story1,obj1] = mklsvm(K,yapp,C,options,verbose);
    timelasso1(i)=toc;
    Kt=mklkernel(xtest,InfoKernel,Weight,options,xapp(posw,:),Sigma);
    ypred=Kt*w+b;
    bc(i)=mean(sign(ypred)==ytest);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [Sigma2,w2,b2,obj2] = compactMKLsvm(K,yapp,C,lambda0,numClusterPos,numClusterNeg,verbose,options);
    timelasso2(i)=toc;
    Kt2=mklkernel(xtest,InfoKernel,Weight,options,xapp,Sigma2);
    ypred2=Kt2*w2+b2;
    bc2(i)=mean(sign(ypred2)==ytest);
end;%



