function play_simple_mkl

%
% play simpleMKL w/ precomputed kernels
%

addpath('../simplemkl1');
addpath('../toollp');

% load global variables
g = exp_setup;



eps         = g.eps;
svmC        = g.svmC;
repetitions = g.repetitions;
nofks       = g.nofks;
nofex       = g.nofex;
cachesize   = g.cachesize;
basis       = g.basis;
ridge_eps   = g.ridge_eps;

% init random number generator
rand('seed',g.rnd_seed);
seed = ceil(rand(repetitions,1)*10000);

fname = sprintf('results/res_simplemkl.mat');
disp(fname);

% result struct:
if(exist(fname)==2)
    load(fname);
else
    smkl.timex    = sparse([]);
    smkl.p_obj    = sparse([]);
    smkl.d_obj    = sparse([]);
    smkl.gap      = sparse([]);
    smkl.globals  = g;
    smkl.p_norm   = 1;
    smkl.solver   = -1;
    smkl.seed     = seed;
    smkl.precompK = 1;
end;

n = nofex;

for k = nofks
    
    for run = 1:repetitions
    
        fprintf('-----------------------------\nrun %d, %d examples, %d kernels\n',run, n, k);
        if(existentry(smkl.timex,k,run)) 
            fprintf('skipping %d kernels, run %d\n',k,run);
            continue;
        end;
        % load data
        [X,y] = mnist_load(n,seed(run));
            
        K = zeros(n,n);
        for j = 1:k
            K(:,:,j) = rbf( X, basis^(j-1) ) + ridge_eps*eye(n);
        end;
    
        %% simple MKL initializations -----------------
        options.algo='svmclass'; 
        % stopping criterion
        options.stopvariation=0;  
        options.stopKKT=0;        
        options.stopdualitygap=1; 
        % choosing the stopping criterion value
        options.seuildiffsigma=1e-2;
        options.seuildiffconstraint=.1;  
        options.seuildualitygap=0.01; %abs(mkl.gap(k,run)/(mkl.d_obj(k,run)));
        % more parameters 
        options.goldensearch_deltmax = 1e-1;
        options.numericalprecision = 1e-8;  
        options.lambdareg = 0;              
        % some algorithms paramaters
        options.firstbasevariable='first';   
        options.nbitermax=100000;         
        options.seuil=0;                     
        options.seuilitermax=10;             
        options.miniter=0;                 
        options.verbosesvm=0;             
        options.efficientkernel=0;           
        %% --------------------------------------------
        
        % run experiment
        tic 
        [beta,w,b,temp,history,obj] = mklsvm(K,y,svmC,options,0);
        smkl.timex(k,run) = toc;
        
        % store result
        smkl.p_obj(k,run) = nan;
        smkl.d_obj(k,run) = obj;
        smkl.gap  (k,run) = nan;
        
        fprintf('t = %ldsec\n', smkl.timex(k,run));
        fprintf('primal_obj = %f, ', smkl.p_obj(k,run));
        fprintf('dual_obj = %f, ', smkl.d_obj(k,run));
        fprintf('gap = %f\n', smkl.gap(k,run));

        save(fname,'smkl');
        
    end;
  
end;

return;
