function play_hessmkl

%
% play hessMKL w/ precomputed kernels
%

addpath('../hessmkl');

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

fname = sprintf('results/res_hessmkl.mat');
disp(fname);

% result struct:
if(exist(fname)==2)
    load(fname);
else
    hmkl.timex    = sparse([]);
    hmkl.p_obj    = sparse([]);
    hmkl.d_obj    = sparse([]);
    hmkl.gap      = sparse([]);
		hmkl.iterations      = sparse([]);
    hmkl.globals  = g;
    hmkl.p_norm   = 1;
    hmkl.solver   = -1;
    hmkl.seed     = seed;
    hmkl.precompK = 1;
end;

% load previous result for stopping criterion
load 'results/res_cplex_1.000norm.mat';


for n = nofex
    
    for run = 1:repetitions
    
        fprintf('-----------------------------\nrun %d, %d examples\n',run, n);
        if(existentry(hmkl.timex,n,run)) 
            fprintf('skipping %d examples, run %d\n',n,run);
            continue;
        end;
        
        % load data
        [X,y] = mnist_load(n,seed(run));
            
        K = zeros(n,n);
        for k = 1:nofks
            K(:,:,k) = rbf( X, basis^(k-1) ) + ridge_eps*eye(n);
        end;
    
        %% simple MKL initializations -----------------
        
        % run experiment
				silp_dgap = abs(mkl.gap(n,run)/(mkl.d_obj(n,run)));
        [beta,t,iterations,p_obj,d_obj,dgap] = kernel_combination(K,y,svmC,silp_dgap);
        hmkl.timex(n,run) = t;
        
        % store result
        hmkl.p_obj(n,run) = d_obj;
				hmkl.d_obj(n,run) = p_obj;
        hmkl.gap(n,run) = dgap;
				hmkl.iterations(n,run) = iterations;
				
        
        fprintf('t = %ldsec\n', hmkl.timex(n,run));
        fprintf('primal_obj = %f, ', hmkl.p_obj(n,run));
        fprintf('dual_obj = %f, ', hmkl.d_obj(n,run));
        fprintf('gap = %f\n', hmkl.gap(n,run));

        save(fname,'hmkl');
        
    end;
  
end;

return;
