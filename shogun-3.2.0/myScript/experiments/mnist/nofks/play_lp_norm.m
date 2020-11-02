function play_lp_norm(p_norm,solver,interleaved)

%
% play MKL w/ precomputed kernels
%
% p_norm >= 1
% solver  = 1 => CPLEX
%         = 2 => Newton/Internal
%

if nargin<3,
	interleaved=1;
end

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

suffix='';
if interleaved==0,
	suffix='_wrapper';
end
switch solver,
case 1,
	fname = sprintf('results/res_cplex_%1.3fnorm%s.mat',p_norm, suffix);
case 2,
	fname = sprintf('results/res_newton_%1.3fnorm%s.mat',p_norm, suffix);
case 3,
    fname = sprintf('results/res_direct_%1.3fnorm%s.mat',p_norm, suffix);
end;
disp(fname);

% result struct:
% result struct:
if(exist(fname)==2)
    load(fname);
else
    mkl.timex    = sparse([]);
    mkl.p_obj    = sparse([]);
    mkl.d_obj    = sparse([]);
    mkl.gap      = sparse([]);
    mkl.globals  = g;
    mkl.p_norm   = p_norm;
    mkl.solver   = solver;
    mkl.seed     = seed;
    mkl.precompK = 1;
end;

n = nofex;

for k = nofks
    
    for run = 1:repetitions
    
        fprintf('-----------------------------\nrun %d, %d examples, %d kernels\n',run, n,k);
        if(existentry(mkl.timex,k,run)) 
            fprintf('skipping %d kernels, run %d\n',k,run);
            continue;
        end;
        % load data
        [X,y] = mnist_load(n,seed(run));
    
        % prepare shogun
        sg('clean_kernel');
        sg('clean_features', 'TRAIN');

        % initialize MKL
        sg('new_classifier', 'MKL_CLASSIFICATION');
	    sg('mkl_use_interleaved_optimization', interleaved); % 0, 1
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
        sg('set_kernel', 'COMBINED', cachesize);
        
        for j = 1:k
            K = rbf( X, basis^(j-1) ) + ridge_eps*eye(n);
            sg('add_kernel', 1, 'CUSTOM', K, 'FULL');
        end;
    
        % run experiment
        tic 
        sg('train_classifier');
        mkl.timex(k,run) = toc;
        
        % store result
        mkl.p_obj(k,run) = sg('compute_mkl_primal_objective');
        mkl.d_obj(k,run) = sg('compute_mkl_dual_objective');
        mkl.gap  (k,run) = mkl.p_obj(k,run) - mkl.d_obj(k,run);
        
        fprintf('t = %ldsec\n',mkl.timex(k,run));
        fprintf('primal_obj = %f, ',mkl.p_obj(k,run));
        fprintf('dual_obj = %f, ',mkl.d_obj(k,run));
        fprintf('gap = %f\n',mkl.gap(k,run));

        save(fname,'mkl');
        
    end;
  
end;

return;
