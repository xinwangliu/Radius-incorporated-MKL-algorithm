function svm = train_sgmklSVM( varargin )

% train_mklSVM - trains a MKL - Suport Vector Machine (SVM) 
%
% Achtung Wichtig!!! Um den Summenkern auf einer Menge S>1 von Gausskernen zu verwenden,
% bitte train_sgmklSVM benutzen!!!
%
% Synopsis:
%     svm = train_mklSVM( properties)
%
% Returns:
%    svm:  a trained SVM struct
%
% Properties:
%        X:    data
%        y:    labels  1 x n    or  2 xn
%        K:   kernel matrix  n x n x k
%        C:   SVM Regularization parameter C
%        mkl_norm:  1<mkl_norm<=inf   (inf==sumKernel (NOT AVERAGE!))
%        X:       {1:n} cell array containing X;   {i}==d_i x n
%        kernel:   'linear' or 'rbf'
%       width:      rbf kernel widths  (e.g. width=logspace(-1,1,5))
%       verbosity: 0, 1 or 2.
%                             If 0: display errors only. (default)
%                            If 1: debug output
%       epsilon:   optimization accuracy (default 1e-3 for SVM 10^(-5) for MKL)
%      'duality_gap': double-check optimization precision in terms of duality gap
%       cachesize: default==100;
%       solver:  'cplex' or 'newton'.
%
%
%==============================================================================
% INIT
%==============================================================================


epsRidge = 1e-10;

properties= propertylist2struct(varargin{:});
properties = set_defaults(properties, ...
                        'y', [], ...
                        'K', [], ...
                        'C',1, ...
                        'verbosity', 0, ...
                        'epsilon', 10^(-3), ...
                        'mkl_norm', [], ...
                        'cachesize', 50, ...
                        'X', [], ...
                        'kernel','linear',...
                        'width', 1,...
                        'duality_gap', inf, ...
                        'solver', 'CPLEX');
struct2workspace(properties);

%  if size(K,3) == 1
%    svm = train_sgSVM(varargin{:});
%    return;
%  end

if (length(width)<=1 && size(K,3)<=1)
  mkl_norm=inf;
end
assert(max(max(size(mkl_norm)))>0);
width = width.^2;

assert(max(max(size(y)))>0);
switch properties.verbosity,
  case 0
    sg('send_command','loglevel ERROR');
    %sg('set_output', '/dev/null');
  case 1
    sg('send_command','loglevel ALL');
  otherwise
    sg('send_command','loglevel ERROR');
end

if mkl_norm==1
  sg('set_solver','CPLEX');
  mkl_epsilon = epsilon*0.01;
else
  if upper(solver(1))=='N' || upper(solver(1))=='I';
    sg('set_solver','INTERNAL');
    mkl_epsilon=epsilon;
  else
    sg('set_solver','CPLEX');
    mkl_epsilon = epsilon*0.01;
  end
end

%prepare shogun
sg('clean_kernel');
sg('clean_features', 'TRAIN');

if size(y,1)==2
  y = sign([y(1,:) + -y(2,:)]);
end
y= sign(y-0.5);
sg('set_labels','TRAIN', y);
sg('set_kernel', 'COMBINED', cachesize);
sg('svm_use_bias', 1);
if mkl_norm<inf
  sg('new_classifier', 'MKL_CLASSIFICATION');
  sg('mkl_parameters', mkl_epsilon, 0, mkl_norm);
  sg('mkl_use_interleaved_optimization', 1); % 0, 1
  sg('set_constraint_generator', 'SVMLIGHT');
else
  sg('new_classifier', 'SVMLIGHT');
  svm.beta = ones(1,max(size(K,3),size(X,3)));
end
sg('svm_epsilon', epsilon);
sg('c', C);

if max(size(K))>0
  K =addRidge(K);
  for k=1:size(K,3)
     sg('add_kernel',1,'CUSTOM', K(:,:,k),'FULL');
  end
else
  assert(max(max(size(X)))>0);
  for i=1:max(size(X,3),length(width))
    sg('add_features','TRAIN', X(:,:,min(i,size(X,3))));
    switch lower(properties.kernel)
      case 'linear'
        sg('send_command', sprintf('add_kernel 1 LINEAR REAL %d %f', cachesize, 1 ));
        svm.kernel{i} = sprintf('LINEAR REAL %d %f', cachesize, 1 );
      case 'gaussian'
        sg('send_command', sprintf('add_kernel 1 GAUSSIAN REAL %d %f', cachesize, width(min(i,length(width))) ));
        svm.kernel{i} = sprintf('GAUSSIAN REAL %d %f', cachesize, width(min(i,length(width)))  );
      otherwise
        error('unknown kernel specified');
      end
    end
end

%sg('send_command', 'set_threshold 0');
sg('train_classifier');
%M=sg('get_kernel_matrix');
[b,alpha]=sg('get_svm') ;
beta = sg('get_subkernel_weights');
svm.beta = beta;
if max(max(size(K)))==0
  svm.suppvec = X(:,alpha(:,2)+1,:);
  alpha(:,2) = [0:size(alpha,1)-1]';
end

svm.alpha=alpha;
svm.b = b;
%obj = abs(sg('get_svm_objective'));
obj=sg('compute_mkl_dual_objective');
svm.duality_gap = 1-obj/sg('compute_mkl_primal_objective');
svm.objective = obj;
svm.mkl_norm=mkl_norm;
  
if svm.duality_gap > properties.duality_gap
  properties.epsilon = properties.epsilon * 0.1;
  varargin = struct2propertylist(properties);
  svm = train_sgmklSVM(varargin{:});
end
  
