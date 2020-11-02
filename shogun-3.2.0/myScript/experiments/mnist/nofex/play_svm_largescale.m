function play_svm_largescale

%
% play SVM = MKL w/ l-inf sum kernels
%

% load global variables
g = exp_setup;
% incorporate additional sample sizes
g.nofex = [g.nofex, g.largescale];
g.cachesize = 3000;

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

fname = sprintf('results/res_svm_ls.mat');
disp(fname);

% result struct:
if(exist(fname)==2)
    load(fname);
else
    svm_ls.timex    = sparse([]);
    svm_ls.p_obj    = sparse([]);
    svm_ls.d_obj    = sparse([]);
    svm_ls.gap      = sparse([]);
    svm_ls.globals  = g;
    svm_ls.p_norm   = inf;
    svm_ls.solver   = -1;
    svm_ls.seed     = seed;
    svm_ls.precompK = 1;
end;


for n = nofex
      
    if (n>=8000)
        repetitions = 1;
    else
        repetitions = g.repetitions;
    end;
    
    for run = 1:repetitions
        
        if(existentry(svm_ls.timex,n,run)) 
            fprintf('skipping %d examples, run %d\n',n,run);
            continue;
        end;
        
        fprintf('-----------------------------\nrun %d, %d examples, %d kernels\n',run, n, nofks);
   
        % load data
        [X,y] = mnist_load(n,seed(run));
    
        % prepare shogun
		%sg('loglevel', 'DEBUG');
        sg('clean_kernel');
        sg('clean_features', 'TRAIN');

        % initialize SVM
        sg('new_svm', 'LIGHT');
        sg('svm_use_bias', 1);
        sg('svm_epsilon', eps);
        sg('c', svmC);
        
        % feed SVM w/ data
        sg('set_labels', 'TRAIN', y');
        sg('set_kernel', 'COMBINED', cachesize);
        
        sg('add_multiple_features',nofks,'TRAIN', X');  
        sg('progress','ON');
        clear X;
        clear y;
        for j = 1:nofks
            sg('add_kernel', 1, 'GAUSSIAN', 'REAL', 0, basis^(j-1));
        end;
               
        % run experiment
        tic 
        sg('train_classifier');
        svm_ls.timex(n,run) = toc;
        
        % store result
        svm_ls.p_obj(n,run) = sg('compute_svm_primal_objective');
        svm_ls.d_obj(n,run) = sg('compute_svm_dual_objective');
        svm_ls.gap  (n,run) = svm_ls.p_obj(n,run) - svm_ls.d_obj(n,run);
        
        fprintf('t = %ldsec\n', svm_ls.timex(n,run));
        fprintf('primal_obj = %f, ', svm_ls.p_obj(n,run));
        fprintf('dual_obj = %f, ', svm_ls.d_obj(n,run));
        fprintf('gap = %f\n', svm_ls.gap(n,run));

        save(fname,'svm_ls');
        
    end;
  
end;

return;
