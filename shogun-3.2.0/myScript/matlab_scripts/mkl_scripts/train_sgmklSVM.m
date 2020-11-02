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
%        elasticnet_lambda:   elastic net parameter,   0<elasticnet_lambda<1
%        kernel:   'linear' or 'rbf'
%       width:      rbf kernel widths  (e.g. width=logspace(-1,1,5))
%       use_bias:   turn bias  on (default) or off
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
                        'degree', 2,...
					    'elasticnet_lambda', [],...
                        'interleaved', true,...
                        'kernel_normalization', 'IDENTITY', ...
                        'linear_kernel_normalization', 'IDENTITY', ...
                        'duality_gap', inf, ...
                        'use_bias',1, ...
                        'svm_solver','LIGHT',...
                        'mkl_solver', 'DIRECT');
struct2workspace(properties);
if exist('normalization')
% RIDGE | VARIANCE | SQRTDIAG | 
  kernel_normalization = upper(normalization);
  linear_kernel_normalization = upper(normalization);
else
  kernel_normalization = upper(kernel_normalization);
  linear_kernel_normalization = upper(linear_kernel_normalization);
end

%  if size(K,3) == 1
%    svm = train_sgSVM(varargin{:});
%    return;
%  end

if length(elasticnet_lambda)==0
  assert(max(max(size(mkl_norm)))>0);
else 
  mkl_solver = 'ELASTICNET';
  mkl_norm = [];
  svm.elasticnet_lambda = elasticnet_lambda;
end
width = 2*width.^2;

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

%for dubugging
%sg('send_command','loglevel DEBUG');
%sg('echo', 'ON');
    
sg('set_solver',upper(mkl_solver));
mkl_epsilon = epsilon*0.01;

if mkl_norm==1 & ~(upper(mkl_solver)=='C' | upper(mkl_solver)=='D')
  mkl_solver = 'DIRECT';
end

%prepare shogun
sg('clean_kernel');
sg('clean_features', 'TRAIN');

if size(y,1)==2
  y = 2*(y(1,:)-0.5);
end
sg('set_labels','TRAIN', y);
if use_bias
  sg('svm_use_bias', 1);
else
  sg('svm_use_bias', 0);
end
sg('set_kernel', 'COMBINED', cachesize);
if (size(K,2)>0 && size(K,3)==1)
  mkl_norm = inf;
end
if length(elasticnet_lambda)>0 || mkl_norm<inf
  sg('new_classifier', 'MKL_CLASSIFICATION');
  if length(elasticnet_lambda)==0
    sg('mkl_parameters', mkl_epsilon, 0, mkl_norm);
  else
    sg('elasticnet_lambda', elasticnet_lambda);
	sg('mkl_parameters', mkl_epsilon, 0);
  end
  sg('mkl_use_interleaved_optimization',interleaved);
  if interleaved
    sg('set_constraint_generator','LIGHT');
  else
    sg('set_constraint_generator',upper(svm_solver));
  end
else
  sg('new_svm', svm_solver);
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
  kernel_str=properties.kernel;
  if ~iscell(kernel_str), 
    aux=kernel_str;
    clear kernel_str;
    kernel_str{1} = aux;
  end
  svm.kernel={};
  svm.kernel_normalization={};
  if size(X,3)>1 && length(kernel_str)==1
    str_aux = kernel_str{1};
    for i=1:size(X,3)
      kernel_str{i}=str_aux;
    end
  end
  if ~(linear_kernel_normalization(1)=='R')
    sg('use_linadd', 1);
  end
  for j=1:length(kernel_str)
    actual_kernel_str=lower(kernel_str{j}(1:4));
    switch lower(actual_kernel_str)
      case 'line'
        sg('add_features','TRAIN', X(:,:,min(j,size(X,3))));
        sg('add_kernel', 1, 'LINEAR', 'REAL', cachesize, 1 );
        svm.kernel{end+1} = sprintf('LINEAR REAL %d %f', cachesize, 1 );
        sg('set_kernel_normalization', kernel_normalization);
        svm.kernel_normalization{end+1}=linear_kernel_normalization;
      case 'gaus'
        sg('use_linadd', 0);
        for i=1:length(width)
          sg('add_features','TRAIN', X(:,:,max(j,size(X,3))));
          sg('add_kernel', 1, 'GAUSSIAN', 'REAL', cachesize, width(i));
          svm.kernel{end+1} = sprintf('GAUSSIAN REAL %d %f', cachesize, width(min(i,length(width))));
          svm.kernel_normalization{end+1}=[];
        end
      case 'poly'
        sg('use_linadd', 0);
        for i=1:length(degree)
          sg('add_features','TRAIN', X(:,:,max(j,size(X,3))));
          inhomogene=false;
          use_normalization=true;
          sg('add_kernel', 1,'POLY', 'REAL', cachesize, degree(i), inhomogene, use_normalization);
          svm.kernel{end+1} = sprintf('POLY REAL %d %f false true', cachesize, width(min(i,length(width))));
          sg('set_kernel_normalization', kernel_normalization);
          svm.kernel_normalization{end+1}=kernel_normalization;
        end
      otherwise
        error('unknown kernel specified');
    end
    sg('init_kernel', 'TRAIN');
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
if mkl_norm==inf
  %obj = abs(sg('get_svm_objective'));
  obj=sg('compute_svm_dual_objective');
  obj2= sg('compute_svm_primal_objective');
else
  obj=sg('compute_mkl_dual_objective');
  obj2= sg('compute_mkl_primal_objective');
end
svm.duality_gap = abs(1-obj/obj2);
svm.objective = obj;
svm.mkl_norm=mkl_norm;
  
if svm.duality_gap > properties.duality_gap
  properties.epsilon = properties.epsilon * 0.1;
  varargin = struct2propertylist(properties);
  set(0,'RecursionLimit',1000);
  if properties.epsilon<10^(-6)
    error('DUALITY GAP IS TOO LARGE!!!  HIGHEST PRECISION IS 10^(-6)');
  end
  svm = train_sgmklSVM(varargin{:});
end

sg('clean_features', 'TRAIN');
sg('clean_kernel');
