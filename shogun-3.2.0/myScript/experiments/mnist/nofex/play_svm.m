function play_svm

%
% play SVM = MKL w/ l-inf sum kernels
%

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

fname = sprintf('results/res_svm.mat');
disp(fname);

% result struct:
if(exist(fname)==2)
    load(fname);
else
    svm.timex    = sparse([]);
    svm.p_obj    = sparse([]);
    svm.d_obj    = sparse([]);
    svm.gap      = sparse([]);
    svm.globals  = g;
    svm.p_norm   = inf;
    svm.solver   = -1;
    svm.seed     = seed;
    svm.precompK = 1;
end;

for n = nofex
    
    for run = 1:repetitions
    
        fprintf('-----------------------------\nrun %d, %d examples\n',run, n);
        
        if(existentry(svm.timex,n,run)) 
            fprintf('skipping %d examples, run %d\n',n,run);
            continue;
        end;
        
        % load data
        [X,y] = mnist_load(n,seed(run));
    
        % prepare shogun
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
        
        K = zeros(n,n);
        for k = 1:nofks
            K = K + rbf( X, basis^(k-1) ) + ridge_eps*eye(n);
        end;
        sg('add_kernel', 1, 'CUSTOM', K, 'FULL');
        
        % run experiment
        tic 
        sg('train_classifier');
        svm.timex(n,run) = toc;
        
        % store result
        svm.p_obj(n,run) = sg('compute_svm_primal_objective');
        svm.d_obj(n,run) = sg('compute_svm_dual_objective');
        svm.gap  (n,run) = svm.p_obj(n,run) - svm.d_obj(n,run);
        
        fprintf('t = %ldsec\n', svm.timex(n,run));
        fprintf('primal_obj = %f, ', svm.p_obj(n,run));
        fprintf('dual_obj = %f, ', svm.d_obj(n,run));
        fprintf('gap = %f\n', svm.gap(n,run));

        save(fname,'svm');
        
    end;
  
end;

return;
