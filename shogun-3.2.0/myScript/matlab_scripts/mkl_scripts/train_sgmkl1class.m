function svm = train_sgmkl1class( varargin )

% train_mklSVM - trains a MKL - Suport Vector Machine (SVM) 
%
% Synopsis:
%     svm = train_sgmkl1class( properties)
%
% Returns:
%    svm:  a trained SVM with fields:
%             alpha: the support vector weights  (sparse format; indexing begins with 0; labels integrated a=y*a)
%             w:      normal vector
%             beta:   the optimal betas
%
% Properties:
%        y:    labels  1 x n    or  2 xn
%        K:   kernel matrix  n x n x k
%        C:   SVM Regularization parameter C
%        mkl_norm:  1<mkl_norm<=inf   (inf==sumKernel (NOT AVERAGE!))
%        linear:       {1:n} cell array;   {i}==d_i x n
%        rbf:            {1:n} cell array;   {i}==d_i x n
%       width:      rbf kernel width  [formula used:  exp((-||.||/width)^2)  ]
%       verbosity: 0, 1 or 2.
%                             If 0: display errors only. (default)
%                            If 1: debug output
%       epsilon:   optimization accuracy (default 1e-5)
%       cachesize: default==250;
%       solver:  'cplex' or 'newton'.
%
%
%==============================================================================
% INIT
%==============================================================================


epsRidge = 1e-10;

properties= propertylist2struct(varargin{:});
properties = set_defaults(properties, ...
                         'K', [], ...
                         'C',1, ...
                         'verbosity', 0, ...
                         'epsilon', 10^(-3), ...
                         'mkl_norm', [], ...
                         'cache_size', 250, ...
                         'linX', [], ...
                         'rbfX', [],...
                         'width', 1,...
                         'solver', 'CPLEX');
struct2workspace(properties);
y = ones(1,size(K,1));
if mkl_norm==inf
  K=mean(K,3);
  mkl_norm = 1;
end
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
else
  if upper(solver(1))=='N';
    sg('set_solver','INTERNAL'); 
  else
    sg('set_solver','CPLEX');
  end
end

%prepare shogun
sg('clean_kernel');
sg('clean_features', 'TRAIN');
sg('set_kernel', 'COMBINED', 0);

% add features/kernels
F=linX;
if max(max(size(F)))>0
  if size(F,1) > 1
    F=F';
  end
  for f=1:size(F,2)
    sg('add_features','TRAIN', F{f});
    sg('add_kernel', 1,  'LINEAR', 'REAL', cache_size);
  end
end
F=rbfX;
if max(max(size(F)))>0
  if size(F,1) > 1
    F=F';
  end
  for f=1:size(F,2)
    sg('add_features','TRAIN', F{f});
    sg('add_kernel', 1, 'GAUSSIAN', 'REAL', cache_size, width(f)^2);
  end
end
if max(size(K))>1
  for k=1:size(K,3)
    sg('add_features','TRAIN',ones(1,size(K,1)));
    sg('add_kernel',1,'CUSTOM', K(:,:,k));
    K =addRidge(K);
    sg('set_custom_kernel', K(:,:,k), 'FULL');
  end
end

if size(y,1)==2
  y = sign([y(1,:) + -y(2,:)]);
end
y= sign(y-0.5);
sg('set_labels','TRAIN', y);
sg('new_svm', 'LIGHT');
sg('use_mkl', 1);
sg('svm_use_bias', false);
sg('mkl_parameters', epsilon, 0, mkl_norm);
sg('svm_epsilon', epsilon);
sg('c', C);
%sg('send_command', 'set_threshold 0');
sg('init_kernel', 'TRAIN');
sg('train_classifier');
M=sg('get_kernel_matrix');
[b,alpha]=sg('get_svm') ;
beta = sg('get_subkernel_weights');
if max(size(linX))>0
  w = [];
  for i=1:size(linX,2)
    w = [w; beta(i)*linX{i}(:,alpha(:,2)+1)*alpha(:,1)];
  end
  svm.w=w;
end

svm.alpha=alpha;
svm.b = b;
svm.beta = beta;
obj = abs(sg('get_svm_objective'));
svm.objective = obj;
