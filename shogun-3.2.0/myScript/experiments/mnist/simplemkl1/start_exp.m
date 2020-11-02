% Example of how to use the mklsvm for  classification
%
%

clear all

addpath('../toollp');
load '../data/data1000.mat'

rand('seed',sum(100*clock));
L = size(K,1);

N = 200;
R = 5;


warning off



verbose=0;

options.algo='svmclass'; % Choice of algorithm in mklsvm can be either
                         % 'svmclass' or 'svmreg'
%------------------------------------------------------
% choosing the stopping criterion
%------------------------------------------------------
options.stopvariation=0; % use variation of weights for stopping criterion 
options.stopKKT=0;       % set to 1 if you use KKTcondition for stopping criterion    
options.stopdualitygap=1; % set to 1 for using duality gap for stopping criterion

%------------------------------------------------------
% choosing the stopping criterion value
%------------------------------------------------------
options.seuildiffsigma=1e-2;        % stopping criterion for weight variation 
options.seuildiffconstraint=0.1;    % stopping criterion for KKT
options.seuildualitygap=0.01;       % stopping criterion for duality gap

%------------------------------------------------------
% Setting some numerical parameters 
%------------------------------------------------------
options.goldensearch_deltmax=1e-3; % initial precision of golden section search
options.numericalprecision=1e-3;   % numerical precision weights below this value
                                   % are set to zero 
options.lambdareg = 1e-10;          % ridge added to kernel matrix 

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
options.efficientkernel=1;         % use efficient storage of kernels 

fprintf('computing.');

for C = [100]
mauc = 0;
    
eps = 10^(-3);

for i=1:R
    fprintf('.');
    ind = randperm(L);
    %train  
    tic
    [beta,w,b,posw,story(i),obj(i)] = mklsvm(K(ind(1:N),ind(1:N),:),y(ind(1:N))',C,options,verbose);
    t(i) = toc;
    
    % classify
    Ktest = Kbeta(K(ind(N+1:end),ind(posw),:),beta',0);
   
    yhat = Ktest*w;
    auc(i) = val_aucs(yhat,y(ind(N+1:end)));
end;



if (mean(auc)>mean(mauc))
    mauc = auc;
    mt   = t;
    mC   = C;
end;

end;


fprintf('done.\n');
fprintf('SimpleMKL (Bach) N = %d, Optimal C = %f\n',N,mC);
fprintf('av-time = %f+-%f seconds\n',mean(mt),std(mt)/sqrt(R));
fprintf('av-auc  = %f+-%f \n',mean(mauc),std(mauc)/sqrt(R));





