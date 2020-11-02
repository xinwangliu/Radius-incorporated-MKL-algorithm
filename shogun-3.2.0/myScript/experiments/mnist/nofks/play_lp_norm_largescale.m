function play_lp_norm_largescale(p_norm,solver);

%
% play MKL w/ feature vectors
%
% p_norm >= 1
% solver  = 1 => CPLEX
%         = 2 => Newton/Internal
%

% check arguments
if ((solver~=1)&&(solver~=2)&&(solver~=3))
    fprintf('please use solver = 1 (CPLEX) or 2 (Newton/Internal) or 3 (Direct)\n');
    return;
end;

if ((p_norm == 1) && (solver==2))
    fprintf('dont use solver 2 (Newton/Internal) for 1 norm mkl\n');
	return
end;

%addpath('~/mkl/src');
%addpath('~/mkl/src/matlab');

% load global variables
g = exp_setup;
% incorporate additional sample sizes
g.nofks = [g.nofks, g.largescale];
g.cachesize = 10000;

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

switch solver,
case 1,
    fname = sprintf('results/res_cplex_%1.3fnorm_ls.mat',p_norm);
case 2,
    fname = sprintf('results/res_newton_%1.3fnorm_ls.mat',p_norm);
case 3,
    fname = sprintf('results/res_direct_%1.3fnorm_ls.mat',p_norm);
end;
disp(fname);


% result struct:
if(exist(fname)==2)
    load(fname);
else
    % result struct:
    mkl_ls.timex    = sparse([]);
    mkl_ls.p_obj    = sparse([]);
    mkl_ls.d_obj    = sparse([]);
    mkl_ls.gap      = sparse([]);
    mkl_ls.globals  = g;
    mkl_ls.p_norm   = p_norm;
    mkl_ls.solver   = solver;
    mkl_ls.seed     = seed;
    mkl_ls.precompK = 0;
end;


n = nofex;

for k = nofks

    if (k>=1000)
        repetitions = 1;
    else
        repetitions = g.repetitions;
    end;
    
    for run = 1:repetitions
    
        if(existentry(mkl_ls.timex,k,run)) 
            fprintf('skipping %d kernels, run %d\n',k,run);
            continue;
        end;
        fprintf('-----------------------------\nrun %d, %d examples, %d kernels\n',run, n, k);
   
        % load data
        [X,y] = mnist_load(n,seed(run));
    
        % prepare shogun
        sg('clean_kernel');
		%sg('loglevel', 'INFO');
		sg('progress', 'ON');
		sg('echo', 'ON');
        sg('clean_features', 'TRAIN');

        % initialize MKL
        sg('new_classifier', 'MKL_CLASSIFICATION');
	    sg('mkl_use_interleaved_optimization', 1); % 0, 1
		sg('set_constraint_generator', 'SVMLIGHT');
        sg('svm_use_bias', 1);
        sg('svm_epsilon', eps);
        sg('c', svmC);
        
        % init MKL
        switch solver,
		case 1,
            % CPLEX
            sg('set_solver','CPLEX');
            mkl_eps = eps*0.01;
		case 2,
            sg('set_solver','NEWTON');
            mkl_eps = eps*0.01;
		case 3,
            sg('set_solver','DIRECT');
            mkl_eps = eps*0.01;
        end;
        sg('mkl_parameters', mkl_eps, 0, p_norm);
         
        % feed SVM w/ data
        sg('set_labels','TRAIN', y');
        sg('set_kernel', 'COMBINED', 0);
              
        sg('add_multiple_features',k,'TRAIN', X'); 
        for j = 1:k  
            sg('add_kernel', 1, 'GAUSSIAN', 'REAL', cachesize/k, basis^(j-1));
        end;

        % run experiment
        tic 
        sg('train_classifier');
        mkl_ls.timex(k,run) = toc;
        
        % store result
        mkl_ls.p_obj(k,run) = sg('compute_mkl_primal_objective');
        mkl_ls.d_obj(k,run) = sg('compute_mkl_dual_objective');
        mkl_ls.gap  (k,run) = mkl_ls.p_obj(k,run) - mkl_ls.d_obj(k,run);
        
        fprintf('t = %ldsec\n',mkl_ls.timex(k,run));
        fprintf('primal_obj = %f, ',mkl_ls.p_obj(k,run));
        fprintf('dual_obj = %f, ',mkl_ls.d_obj(k,run));
        fprintf('gap = %f\n',mkl_ls.gap(k,run));

        save(fname,'mkl_ls');
        
    end;
  
end;

return;
