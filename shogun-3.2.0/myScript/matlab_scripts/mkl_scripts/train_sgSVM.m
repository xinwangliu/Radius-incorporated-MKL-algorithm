function sgSVM = train_sgSVM( varargin);
% train_sgSVM - trains a Suport Vector Machine (SVM) using the genefinder program
%
%  IMPORTANT!!! THIS IS AN OBSOLETE FUNCTION USE train/apply_sgmklSVM !!!
% This function is not longer updated!!
%
% Synopsis:
%   sgSVM = train_sgSVM(data, labels, properties)
%   sgSVM = train_sgSVM(data, labels, 'Property', Value, ...)
%   
% Arguments:
%   data:	[d,n]	real valued training data containing n examples of d dims
%	labels:	[2,n]	logical labels (2-class problems only!) or [1 d]
%                       vector with +1/-1 entries
%
% Returns:
%   sgSVM:	a trained SVM with fields:
%             'SV':     the support vectors
%             'alphas': the alphas for each support vector
%             'w':      normal vector
%             'b':      the bias
%             'kernel': the kernel with parameters as string
%
% Properties:
%	'C': scalar or [1 2] matrix. Regularization parameter C
%                       If C is a [1 2] matrix, class wise regularization
%                       will be used, with C(1) the parameter for class 1
%                       (the "-1 class") and C(2) for class 2 (the "+1 class")
%
%       'C_weight'      weight by which the regularization parameter C will 
%                       be multiplied for the smaller class (default 1). 
%                       Only used when length(C) == 1
%
% 	'kernel': the name of the kernel to be used
%			linear: K(u,v) = u'*v;
%			poly: K(u,v) = (u'*v)^d when homogeneous
%				  K(u,v) = (u'*v+1)^d when inhomogeneous
%			gaussian: K(u,v) =  exp(-||u-v||^2/width)
%
%	'implementation': Which SVM implementation to use:
%			light: SVMLight (see http://svmlight.joachims.org) (default)
%			libsvm: libSVM (see http://www.csie.ntu.edu.tw/~cjlin/libsvm )
%     svmocas: Voitecs SVM
%
%       'verbosity'     0, 1 or 2. 
%                       If 0: display errors only.
%                       If 1: display errors and warnings (Default)
%                       If 2: debug output
%	'width':	gaussian kernel width
%	'degree:'	degree of the polynomial
%	'epsilon':  optimization accuracy (default 1e-3)
%    'duality_gap': double-check optimization precision in terms of duality gap
%	'cachesize': size of the kernel cache in megabytes (default 100MB)
%	'poly_inhomogene:'	if set to 0 homogene polynomial kernel, set to 1 for
%						inhomogene kernel
% 'K':  custom training kernel matrix
%
% Description:
%   The function train_SVM trains a Support Vector Machine, i.e. it
%	determines the support vectors x_i, lagrange multipliers alpha_i and
%	the bias term b of the following function
%		f(x)=sum_{i=0}^l alpha_i K(x_i,x) + b.
%	The training algorithm maximizes the margin between the classes and
%	solves a quadratic program (see References for more details).
%
% Examples:
%   x = [randn(10,50)-0.5, randn(10,50)];
%	y = [repmat([1;0], 1, 50) repmat([0;1], 1,50)];
%
%
%   s = train_sgSVM(x,y, 'C',1, 'kernel','linear', 'implementation','light');
%			trains a standard linear SVM using SVMlight. not very useful.
%
%   s = train_sgSVM(x,y, 'C',5, 'kernel','poly', 'degree', 3, 'svm','libsvm');
%			trains an inhomogenious polynomial svm of degree 3 using LibSVM
%
%	c=train_sgSVM(x, y, 'width', 10, 'cachesize', 100, 'C', [1 2] , 'kernel', 'gaussian', 'implementation', 'light');
%	out=apply_sgSVM(c, x);
%			trains a rbf svm with sigma^2=10 different Cs for each class
%			using SVMlight
%
% References:
%   Vapnik:The Nature of Statistical Learning Theory. Springer, 1995.
%   Burges: A Tutorial on Support Vector Machines for Pattern Recognition, 1998
%
% See also: apply_sgSVM

% Soeren Sonnenburg (2005)
% modified S. Mika (2005)
% based on work from die Guido Dornhege (2004)
% $Id: train_sgSVM.m,v 1.13 2005/06/10 10:46:06 neuro_toolbox Exp $


%error(nargchk(3, 100, nargin));
warning off;

%% bb: setting up properties in the nice way
properties= propertylist2struct(varargin{:});
properties= set_defaults(properties, ...
                         'X', [], ...
                         'y', [], ...
                         'C',1, ...
                         'kernel','linear', ...
                         'implementation', 'libsvm', ...
                         'width', 1, ...
                         'degree', 1, ...
                         'cachesize', 100, ...
                         'verbosity', 0, ...
                         'poly_inhomogene', 0, ...
                         'epsilon', 1e-3, ...
                         'duality_gap', inf, ...
                         'C_weight', 1, ...
                         'trainloss', 0, ...
                         'use_bias', 1, ...
                         'K', [] );

struct2workspace(properties);

assert(max(max(size(y)))>0);
if size(y,1)>size(y,2), y=y'; end
if size(y,1)==1,
  % y as +1/-1 vector. Do not do any error checking, just make sure
  % it is really only +1/-1
  y = sign(y);
elseif size(y,1)==2,
  y = [-1 1]*y;
else
  error('sgSVM only works for binary classification tasks');
end

%Test von SVM OCAS
%implementation = 'svmocas';

if exist('Ktr')&& max(max(size(K)))==0
  K = Ktr;
end
Ktr = [];

if size(K,2)>1 && size(X,1)==0
  X = ones(1,size(K,2));
end

switch verbosity,
  case 0
    sg('send_command','loglevel ERROR');
    %sg('set_output', '/dev/null');
  case 1
    sg('send_command','loglevel ALL');
end

if use_bias
  sg('svm_use_bias', 1);
else
  sg('svm_use_bias', 0);
end

switch lower(properties.kernel),
  case 'linear',
    sgSVM.kernel_string = ['LINEAR REAL ', int2str(cachesize), ' 1.0'];
  case 'poly',
    sgSVM.kernel_string = ['POLY REAL ', int2str(cachesize), ' ', int2str(degree), ' ', num2str(poly_inhomogene)];
  case 'gaussian',
    sgSVM.kernel_string = ['GAUSSIAN REAL ',  int2str(cachesize), ' ', num2str(width)];
  otherwise
    if size(K,1)==0
      error('unknown kernel specified');
    end
    sgSVM.kernel_string = 'SWITCH_TO_CUSTOM_KERNEL';
end

switch lower(implementation)
  case 'light',
    implementation = 'LIGHT'
    sg('set_features','TRAIN', full(X));
  case 'libsvm',
    implementation = 'LIBSVM';
    sg('set_features','TRAIN', full(X));
  case 'svmocas',
    implementation = 'SVMOCAS';
    if size(K>0)
      X = empKM(K);
      K = [];
    end
    sg('set_features','TRAIN', sparse( [ X ; ones(1,size(X,2)) ] ));
  otherwise
    implementation = 'LIGHT';
    sg('set_features','TRAIN', full(X));
end

% $$$ sg('send_command','set_output /dev/null');
  

sg('set_labels','TRAIN', y);

if implementation(4) == 'o'
  sg('svm_use_bias', 0);
  sg('svm_bufsize', 1000);
  sg('new_classifier', 'SVMOCAS');
else
  if size(K,1)==0
    sg('send_command', sprintf('set_kernel %s',sgSVM.kernel_string));
  else
    sg('set_kernel', 'CUSTOM', addRidge(K), 'FULL');
end
  sg('send_command', 'init_kernel TRAIN');
  sg('send_command', sprintf('new_svm %s', implementation));
end

% Reweightning of class weights
if (C_weight ~= 1) & (length(C) == 1)
  % C_weight-times higher value of parameter C for smaller class
  n_size = sum(y<0);
  p_size = sum(y>0);
  C = [1,1] * C;
  if n_size<p_size
    C(1) = C(1) * C_weight;
  else
    C(2) = C(2) * C_weight;
  end
end
sg('send_command', sprintf('c%s', sprintf(' %f', C)));
sg('send_command', sprintf('svm_epsilon %s', epsilon));
sg('send_command', 'svm_train');

[sgSVM.b, sgSVM.alphas]=sg('get_svm');

if size(K,2)==0
  switch lower(implementation)
    case 'svmocas',
      [foo, wb]=sg('get_svm');
      sgSVM.w = wb (1:end-1,:);
      sgSVM.b =  wb(end,1);
    otherwise
      [sgSVM.b, sgSVM.alphas]=sg('get_svm');
      sgSVM.suppvec = X(:,sgSVM.alphas(:,2)+1);
      sgSVM.alphas(:,2) = [0:size(sgSVM.alphas,1)-1]';
  end
  sgSVM.w=X(:,sgSVM.alphas(:,2)+1)*sgSVM.alphas(:,1);
  if max(size(K)>0)
    sgSVM.objVal = sum(sgSVM.w.^2)+C*sum((1-y.*(sgSVM.w'*X+sgSVM.b)).^2);
    sgSVM.alphas = pinv(X)*sgSVM.w;
    sgSVM.alphas(:,2) = 0:size(sgSVM.alphas,1)-1;
    sgSVM.suppvec = X(:,sgSVM.alphas(:,2)+1);
  end
  if trainloss
    out = apply_sgSVM(sgSVM, X);
    sgSVM.trainloss = mean( y~= sign(out));
  end
else
   % Achtung in den Alpha steck schon die y drin: alpha_y_i = y_i * alpha_i !!
   [sgSVM.b, sgSVM.alphas]=sg('get_svm');
   alphaY = sgSVM.alphas(:,1);
   svidx = sgSVM.alphas(:,2)+1;
   sgSVM.objVal = sum(abs(alphaY)) - 0.5*alphaY' * K(svidx,svidx) * alphaY;
   if trainloss
     out = apply_sgSVM(sgSVM, [], 'Kts', K);
   end
end

if trainloss
  %y = abs([y+1; y-1])/2;
  %sgSVM.trainErr_0_1 = 1-mean(loss_0_1(y,out));
end

sgSVM.epsilon=epsilon;
sgSVM.duality_gap = 1-sg('compute_svm_dual_objective')/sg('compute_svm_primal_objective');

if sgSVM.duality_gap > duality_gap
  properties.epsilon = epsilon * 0.1;
  varargin = struct2propertylist(properties);
  sgSVM = train_sgSVM(varargin{:});
end