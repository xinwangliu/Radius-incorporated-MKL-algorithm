function result = exp_mkl_toy_gap

%--------------------------------------
%  Parameters
%--------------------------------------

plot_gap = 1;

gaps=logspace(-3,0,7);
ntrs = 50;
nte = 1000;
nval = 500;
reps = 100;
Cs = logspace(-3,1,9);
N = 50;
d =1;
sep = 4;
verbosity = 0;
epsilon = 10^(-3);
ridge = 10^(-10);
mkl_norms = [1,  2, inf];
duality_gap = 0.01;

%--------------------------------------
% Assertions & Cases
%--------------------------------------
tic;
if plot_gap
  assert(max(size(ntrs))==1);
end
n = max(ntrs) + nval + nte;
assert(mod(n,2)==0);
assert(all(mod(ntrs,2)==0));
assert(all(gaps>0));

%--------------------------------------
% Main Loops
%--------------------------------------

tic; t=0;
for a=1:length(gaps)
  gap = gaps(a);

  %-------------------
  % Prepare mu's
  %-------------------
  block_mu = fvnorm(exp(-[1:gap:1+(N/d-1)*gap]))';
  block_mu = block_mu/sum(abs(block_mu));
  mu = fvnorm(reshape(repmat(block_mu',d,1),N,1));
  mu = sep * mu / norm(mu);
  block_mus(:,a) = block_mu;
  mus(:,a) = mu;

  for r =1:reps
    %-------------------
    % generate toy data
    %-------------------
    selector = repmat([-1,1],1,max(ntrs)/2);
    ytr=sign(selector-0.5);
    yte=repmat( [1,-1], 1, (nval+nte)/2);
    y = [ ytr, yte ];
    X =[ randn(N,max(ntrs))+ (ones(N,1)*sign(selector-0.5)).*(0.5*mu*ones(1,max(ntrs))) , randn(N,nval+nte) + repmat( [0.5*mu,-0.5*mu] ,1 ,(nval+nte)/2 ) ];
    
    %-------------------
    % calculate kernels for gap/ntr experiment
    %-------------------
    X=reshape(X', [1 fliplr(size(X))]);
    for ni = 1:length(ntrs)
      %-------------------
      % divide
      %-------------------
      ntr = ntrs(ni);
      divtr = [1:ntr];
      divval = [max(ntrs)+1:max(ntrs)+nval];
      divte = [max(ntrs)+nval+1:max(ntrs)+nval+nte];

      for m = 1:length(mkl_norms)
        mkl_norm = mkl_norms(m);
        for c = 1:length(Cs)
          %-------------------
          % train
          %-------------------
          C = Cs(c);
          t=t+1;
          svm = train_sgmklSVM( 'X', X(1,divtr,:), 'y', y(1,divtr), 'C',C,'verbosity',verbosity, 'mkl_norm',mkl_norm,'epsilon',epsilon,'normalization','VARIANCE','kernel','linear','mkl_solver','DIRECT');
          duality_gaps(a,m,c,r,ni) = svm.duality_gap;
          outval = apply_sgmklSVM( 'svm', svm, 'X', X(1,divval,:),'verbosity',verbosity, 'mkl_norm', mkl_norm,'epsilon',epsilon);
          outte = apply_sgmklSVM( 'svm', svm, 'X', X(1,divte,:),'verbosity',verbosity, 'mkl_norm', mkl_norm, 'epsilon',epsilon');
          errs_val(a,ni,m,c,r) = mean( y(1,divval) ~= sign(outval) );
          errs(a,ni,m,c,r) = mean( y(1,divte) ~= sign(outte) );
          delta_beta(a,ni,m,c,r) = norm( svm.beta'/norm(svm.beta) - block_mu/norm(block_mu) );
          print_progress(t, length(gaps)*length(ntrs)*length(mkl_norms)*length(Cs)*reps);
        end
        
      end
      
    end
    
  end

end