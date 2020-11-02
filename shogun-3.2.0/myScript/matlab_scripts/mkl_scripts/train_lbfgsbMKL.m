function svm = train_lbfgsbMKL( varargin )

% train_bfgsSVM - trains a MKL - elasticnet/lp-norm - Suport Vector Machine (SVM) with L-BFGS
%
%
% Synopsis:
%     svm = train_bfgsSVM( properties)
%
% Returns:
%    svm:  a trained SVM struct
%
% Properties:
%        y:    labels  1 x n    or  2 xn
%        K:   kernel matrix  n x n x k
%        C:   SVM Regularization parameter C
%        lambda:  elasticnet parameter
%        mkl_norm:  1<mkl_norm<=inf   (inf==sumKernel (NOT AVERAGE!))
%
%
%
%  We solve the formulation:  C\sum_i\loss_i + 1/2*||w||_{2,p}^2 + \lambda/2*||w||_2^2
%
% Author: Marius Kloft    (April, 2010)
%

properties= propertylist2struct(varargin{:});
properties = set_defaults(properties, ...
                        'y', [], ...
                        'K', [], ...
                        'C',1, ...
						'lambda',0, ...
                        'mkl_norm', []);
struct2workspace(properties);


assert(max(max(size(y)))>0);
if size(y,1)==2
  y = 2*(y(1,:)-0.5);
end


if lambda>0
  if exist('mkl_norm')
    %assert(mkl_norm<1.05);
  else 
    mkl_norm=1.02;
  end
end
assert(max(max(size(mkl_norm)))>0);
if mkl_norm == 1 
  mkl_norm=1.02;
end

ridge=10^(-10);
K =addRidge(K, mean(diag(mean(K,3)))*10^(-10));


if ~(mkl_norm<inf && mkl_norm>1)
  error('mkl_norm should be 1<p<inf');
else
  %here goes the mkl training!
  alpha0 = 0.1*C*ones(size(K,1),1);
  gamma0 = zeros(size(K,1),1);
  lb = {zeros(size(K,1),1)  -inf*ones(size(K,1),1)}; 
  ub = {C*ones(size(K,1),1)  inf*ones(size(K,1),1)};
  [alpha gamma] = lbfgsb({alpha0 gamma0},lb,ub,'computeObjectiveLbfgsbMKL','computeGradientLbfgsbMKL',{K,y',lambda,mkl_norm},'nocallback','maxiter',3000,'m',4,'factr',1e-12,'pgtol',1e-5);
  %[alpha gamma] = lbfgsb({alpha0 gamma0},lb,ub,'computeObjectiveLbfgsbMKL','computeGradientLbfgsbMKL',{K,y',lambda,mkl_norm},'genericcallback','maxiter',5000,'m',4,'factr',1e-12,'pgtol',1e-4);
  svm.alpha = alpha'.*y;
  svm.gamma = gamma;
  svm.primal_objective = computeObjectiveLbfgsbMKL(alpha,gamma,{K,y',lambda,mkl_norm});
  
  assert(~isnan(svm.primal_objective));
  
  quadterm = zeros(1,size(K,3));
  for m=1:size(K,3)
	quadterm(m) = ((svm.alpha')'*K(:,:,m)*(svm.alpha'))^(0.5);
  end
  if (lambda>0)
    factor = 1 - lambda^(mkl_norm-1);
    %svm.beta_alt = max(0,(lambda*((((lambda+size(K,3))*quadterm ) / sum(quadterm) -1).^(-1)) + lambda).^(-1));
	  beta = 1/lambda*max(0, 1-sum(quadterm)./((lambda+size(K,3))*quadterm));
	%beta = 1/lambda*max(0, 1-sum(quadterm)./((lambda+sum(beta.^(mkl_norm-1)))*quadterm));
	beta = 1/lambda*max(0, 1-sum(quadterm)./((lambda+sum(beta==0))*quadterm));
	%beta = 1/lambda*max(0, 1-beta.^(mkl_norm-1)*sum(quadterm)^(2-mkl_norm)./((lambda+sum(beta.^(mkl_norm-1)))^(2-mkl_norm)*quadterm.^(2-mkl_norm))).*logical(beta~=0);
	%beta = max(0,lambda^(-1) * (1-sum(quadterm.^(0.5))./((size(K,3)+lambda)*quadterm.^(0.5))));
	if sum(beta)==0, 
	  beta = 1/lambda*max(0, 1-sum(quadterm)./((lambda+size(K,3))*quadterm));
	end
	alpha=alpha(:); gamma=gamma(:); y=y(:);
	for i=1:size(K,3)
	 quadterm(i) = (alpha.*y-gamma)'*K(:,:,i)*(alpha.*y-gamma);
	end
	for i=1:size(K,3)
	 kernel_isActive(i) = (abs(1-quadterm(i)/max(quadterm))<0.1);
	end
	M_active = sum(kernel_isActive);
	beta_tilde = 1/lambda*max(0, ((lambda+M_active)*quadterm)./sum(quadterm.*kernel_isActive)-1);
	beta = beta_tilde;
	svm.beta = (beta_tilde.^(-1)+lambda).^(-1);
  else
    svm.beta = quadterm.^((2-mkl_norm)/(mkl_norm-1));
  end
end
svm.beta =svm.beta /sum(svm.beta);


% if ~(mkl_norm<inf && mkl_norm>1)
  % error('mkl_norm should be 1<p<inf');
% else
  % %here goes the mkl training!
  % alpha0 = zeros(size(K,1),1);
  % beta0 = zeros(size(K,1),1);
  % if lambda >0.0000000001
    % lb = {zeros(size(K,1),1) , -inf*ones(size(K,1),1)};
    % ub = {C*ones(size(K,1),1) , inf*ones(size(K,1),1)};
  % else
    % lb = {zeros(size(K,1),1) , zeros(size(K,1),1)};
    % ub = {C*ones(size(K,1),1) , zeros(size(K,1),1)};
  % end
  % [alpha, beta]= lbfgsb({alpha0, beta0},lb,ub,'computeObjectiveLbfgsbMKL','computeGradientLbfgsbMKL',{K,y,lambda,mkl_norm});
       % % ...,'genericcallback','maxiter',80,'m',4,'factr',1e-12,'pgtol',1e-5);
  % svm.alpha = alpha,
  % svm.beta = beta;
% end




%svm.kernel={};
%svm.kernel_normalization={};
 

%svm.beta = beta;
%  alpha(:,2) = [0:size(alpha,1)-1]';
% obj=sg('compute_svm_dual_objective');
%  obj2= sg('compute_svm_primal_objective');
%svm.duality_gap = abs(1-obj/obj2);
svm.mkl_norm=mkl_norm;
svm.b=0;  


