function sgSVM = train_sg1classSVM(data, labels, varargin);
% train_sgSVM - trains a Suport Vector Machine (SVM) using the genefinder program
%
% Synopsis:
%   sgSVM = train_sg1classSVM(data, labels, properties)
%   sgSVM = train_sg1classSVM(data, labels, 'Property', Value, ...)
%   
% Arguments:
%   data:	[d,n]	real valued training data containing n examples of d dims
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
%	'epsilon':  optimization accuracy (default 1e-5)
%	'cachesize': size of the kernel cache in megabytes (default 100MB)
%	'poly_inhomogene:'	if set to 0 homogene polynomial kernel, set to 1 for
%						inhomogene kernel
% 'Ktr':  custom training kernel matrix
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
                         'C',1, ...
                         'kernel','linear', ...
                         'width', 10, ...
                         'degree', 1, ...
                         'cachesize', 250, ...
                         'verbosity', 0, ...
                         'poly_inhomogene', 0, ...
                         'epsilon', 1e-3, ...
                         'C_weight', 1, ...
                         'trainloss', 0, ...
                         'Ktr', [], ...
                         'K', [] );

%Test von SVM OCAS
%properties.implementation = 'svmocas';

if size(properties.Ktr,1)==0
  properties.Ktr = properties.K;
end
properties.K = [];

if size(properties.Ktr,2)>1 && size(data,1)==0
  data = ones(1,size(properties.Ktr,2));
end

labels = ones(1,size(properties.Ktr,1));

switch properties.verbosity,
  case 0
    sg('send_command','loglevel ERROR');
    %sg('set_output', '/dev/null');
  case 1
    sg('send_command','loglevel ALL');
end


switch lower(properties.kernel),
  case 'linear',
    kernel = ['LINEAR REAL ', int2str(properties.cachesize), ' 1.0'];
  case 'poly',
    kernel = ['POLY REAL ', int2str(properties.cachesize), ' ', int2str(properties.degree), ' ', num2str(properties.poly_inhomogene)];
  case 'gaussian',
    kernel = ['GAUSSIAN REAL ',  int2str(properties.cachesize), ' ', num2str(properties.width)];
  otherwise
    error('unknown kernel specified');
end

sgSVM.kernel_string=kernel;

implementation = 'LIGHT';
sg('set_features','TRAIN', full(data));

sg('set_labels','TRAIN', labels);
sg('svm_use_bias', 0);

if size(properties.Ktr,1)==0
  sg('send_command', sprintf('set_kernel %s',sgSVM.kernel_string));
else
  sg('set_kernel', 'CUSTOM');
  sg('set_custom_kernel', addRidge(properties.Ktr), 'FULL');
end
sg('send_command', 'init_kernel TRAIN');
sg('send_command', sprintf('new_svm %s', implementation));


% Reweightning of class weights
if (properties.C_weight ~= 1) & (length(properties.C) == 1)
  % C_weight-times higher value of parameter C for smaller class
  n_size = sum(labels<0);
  p_size = sum(labels>0);
  properties.C = [1,1] * properties.C;
  if n_size<p_size
    properties.C(1) = properties.C(1) * properties.C_weight;
  else
    properties.C(2) = properties.C(2) * properties.C_weight;
  end
end
sg('send_command', sprintf('c%s', sprintf(' %f', properties.C)));
sg('send_command', sprintf('svm_epsilon %s', properties.epsilon));
sg('send_command', 'svm_train');

[sgSVM.b, sgSVM.alphas]=sg('get_svm');

if size(properties.Ktr,2)==0
  switch lower(properties.implementation)
    case 'svmocas',
      [foo, wb]=sg('get_svm');
      sgSVM.w = wb (1:end-1,:);
      sgSVM.b =  wb(end,1);
    otherwise
      [sgSVM.b, sgSVM.alphas]=sg('get_svm');
  end
  sgSVM.w=data(:,sgSVM.alphas(:,2)+1)*sgSVM.alphas(:,1);
  if max(size(properties.Ktr)>0)
    sgSVM.objVal = sum(sgSVM.w.^2)+properties.C*sum((1-labels.*(sgSVM.w'*data+sgSVM.b)).^2);
    sgSVM.alphas = pinv(data)*sgSVM.w;
    sgSVM.alphas(:,2) = 0:size(sgSVM.alphas,1)-1;
    sgSVM.suppvec = data(:,sgSVM.alphas(:,2)+1);
  end
  if properties.trainloss
    out = apply_sgSVM(sgSVM, data);
  end
else
   % Achtung in den Alpha steck schon die labels drin: alpha_y_i = y_i * alpha_i !!
   [sgSVM.b, sgSVM.alphas]=sg('get_svm');
   alphaY = sgSVM.alphas(:,1);
   svidx = sgSVM.alphas(:,2)+1;
   sgSVM.objVal = sum(abs(alphaY)) - 0.5*alphaY' * properties.Ktr(svidx,svidx) * alphaY;
   if properties.trainloss
     out = apply_sgSVM(sgSVM, [], 'Kts', properties.Ktr);
   end
end

if properties.trainloss
  labels = abs([labels+1; labels-1])/2;
  sgSVM.trainErr_0_1 = 1-mean(loss_0_1(labels,out));
end

sgSVM.epsilon=properties.epsilon;

