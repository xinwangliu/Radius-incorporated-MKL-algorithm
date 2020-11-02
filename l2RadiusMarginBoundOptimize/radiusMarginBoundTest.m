clear;
clc;
% Example of how to use the mklsvm for  classification
%
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
options.numericalprecision=1e-8;   % numerical precision weights below this value
                                   % are set to zero 
options.lambdareg = 1e-8;          % ridge added to kernel matrix 

%------------------------------------------------------
% some algorithms paramaters
%------------------------------------------------------
options.firstbasevariable='first'; % tie breaking method for choosing the base 
                                   % variable in the reduced gradient method 
options.nbitermax=500;             % maximal number of iteration  
options.seuil=0;                   % forcing to zero weights lower than this 
options.seuilitermax=5;           % value, for iterations lower than this one 

options.miniter=0;                 % minimal number of iterations 
options.verbosesvm=0;              % verbosity of inner svm algorithm 
options.efficientkernel=0;         % use efficient storage of kernels 


%------------------------------------------------------------------------
%                   Building the kernels parameters
%------------------------------------------------------------------------
kernelt={ 'gaussian' 'poly'};
kerneloptionvect={ [0.5 1 2 5 7 10 12 15 17 20] [1:10]}; %%[0.5 1 2 5 7 10 12 15 17 20]  
variablevec={ 'all' 'all'};  %'single'
ratio  = 0.6;
nbiter = 30;
classcode=[1 -1];
lambda0 = options.lambdareg;

cSet = 2.^[-5:2:15];

path = 'F:\XINWAN_WORK\MKLToolBox\mkLearningWithRadius\dataset\';
dataName = 'wpbc';
load([path dataName, '.mat']);

[nbdata,dim]=size(trnData);
nbtrain=floor(nbdata*ratio);
for i=1: nbiter
    i
    [xapp,yapp,xtest,ytest,indice]=CreateDataAppTest(trnData, trnLabel, nbtrain,classcode);
    save(['F:\XINWAN_WORK\MKLToolBox\mkLearningWithRadius\10GroupData\wpbc\',dataName,'_generalization_',num2str(i),'.mat'],'xapp', 'yapp', 'xtest', 'ytest','indice');
%     [xapp,xtest]=normalizemeanstd(xapp,xtest);
%     [kernel,kerneloptionvec,variableveccell]=CreateKernelListWithVariable(variablevec,dim,kernelt,kerneloptionvect);
%     [Weight,InfoKernel]=UnitTraceNormalization(xapp,kernel,kerneloptionvec,variableveccell);
%     
%     K=mklkernel(xapp,InfoKernel,Weight,options);
% %     
% % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     C = 500;
%     disp('simpleMKL');
%      tic
%     [beta,w,b,posw,story(i),obj(i)] = mklsvm(K,yapp,C,options,verbose);
%     timelasso(i)=toc
% 
%     Kt=mklkernel(xtest,InfoKernel,Weight,options,xapp(posw,:),beta);
%     ypred=Kt*w+b;
% 
%     bc1(i)=mean(sign(ypred)==ytest)
% 
% % % %     
% % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp('MKL Radius Margin Bound--Tuning C by optimizing')
%     tic
%     [beta2,w2,b2,posw2,story2(i),obj2(i)] = radiusMarginBoundVersion2(K,yapp,options,verbose);
%     timelasso(i)=toc
% 
%     tt2(i,:) = [1/beta2(1) sum(beta2(2:end))];
%     Kt2=mklkernel(xtest,InfoKernel,Weight,options,xapp(posw2,:),beta2(2:end));
%     ypred2=Kt2*w2+b2;
% 
%     bc2(i)=mean(sign(ypred2)==ytest)
% %     
% % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp('optimal Margin -- Tuning C by optimization');
%     tic
%     [beta3,w3,b3,posw3,story3(i),obj3(i)] = l2MarginMKL(K,yapp,options,verbose);
%     timelasso(i)=toc
% 
%     tt3(i,:) = [1/beta3(1) sum(beta3(2:end))];
%     Kt3=mklkernel(xtest,InfoKernel,Weight,options,xapp(posw3,:),beta3(2:end));
%     ypred3=Kt3*w3+b3;
% 
%     bc3(i)=mean(sign(ypred3)==ytest)
%     
% %     C1=100;
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% [Sigma,Alpsup,w0,pos,history,obj,status] = l2RadiusMarginBoundCV(K,yapp,C,option,verbose)
%     disp('MKL Radius Margin Bound -- Tuning C by CV');
%      tic
%     [beta4,w4,b4,posw4,story4(i),obj4(i)] = l2RadiusMarginBoundCV(K,yapp,C1,options,verbose);
%     timelasso(i)=toc
%     
%     Kt4=mklkernel(xtest,InfoKernel,Weight,options,xapp(posw4,:),beta4);
%     ypred4=Kt4*w4+b4;
%     
%     bc4(i)=mean(sign(ypred4)==ytest)
%     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%    %% [Sigma,Alpsup,w0,pos,history,obj,status] = l2optimalMarginCV(K,yapp,C,option,verbose);
%     
%    disp('optimal Margin -- Tuning C by CV');
%    tic
%    for iter = 1:length(cSet)
%        [beta5,w5,b5,posw5,story5,obj5] = l2optimalMarginCV(K,yapp,cSet(iter),options,verbose);
%        timelasso(i)=toc
%        
%        Kt5=mklkernel(xtest,InfoKernel,Weight,options,xapp(posw5,:),beta5);
%        ypred5=Kt5*w5+b5;
%        
%        bc5(iter)=mean(sign(ypred5)==ytest)
%        
%        margin(iter) = obj5;
%    end
end



