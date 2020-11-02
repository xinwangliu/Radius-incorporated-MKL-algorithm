function play_heshmkl

%
% play heshmkl w/ precomputed kernels
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

n = nofex;

for k = nofks
    
    for run = 1:repetitions
    
        fprintf('-----------------------------\nrun %d, %d examples, %d kernels\n',run, n, k);
        if(existentry(hmkl.timex,k,run)) 
            fprintf('skipping %d kernels, run %d\n',k,run);
            continue;
        end;
        % load data
        [X,y] = mnist_load(n,seed(run));
            
        K = zeros(n,n);
        for j = 1:k
            K(:,:,j) = rbf( X, basis^(j-1) ) + ridge_eps*eye(n);
						%K(:,:,j) = center(scale(K(:,:,j)));
        end;
    

        % run experiment
				%silp_dgap = -1; % abs(mkl.gap(n,run)/(mkl.d_obj(n,run)));
				silp_dgap = 0.01;
				[beta,t,iterations,p_obj,d_obj,dgap] = kernel_combination(K,y,svmC,silp_dgap);
        hmkl.timex(k,run) = t;

        
        % store result
        hmkl.p_obj(k,run) = d_obj;
				hmkl.d_obj(k,run) = p_obj;
        hmkl.gap(k,run) = dgap;
				hmkl.iterations(k,run) = iterations;
        
        fprintf('t = %ldsec\n', hmkl.timex(k,run));
        fprintf('primal_obj = %f, ', hmkl.p_obj(k,run));
        fprintf('dual_obj = %f, ', hmkl.d_obj(k,run));
        fprintf('gap = %f\n', hmkl.gap(k,run));

        save(fname,'hmkl');
        
    end;
  
end;

return;
