function [D,t,iteration,p_obj,d_obj,dgap] = kernel_combination(K_,Y,C,max_dgap)
%
% D = KERNEL_COMBINATION(K,Y)
% Learn a combination of kernel given by the 3D matrix K
% Postive weights are returned in D [note that length(D) = size(K,3)]
%
% The kernels need to be centered and normalized
%
% C is the constant penalizing linearly the training errors.
% If C is not given, but one of the kernels is the identity, this will 
% autmatically learn the right C with an L2 penalization of the slacks
%
% Requires the SVQP quadratic solver of Leon Bottou.
% Can also work with the quadprog of Matlab but it is slower.
% 
% Details of the algorithm available at: 
% http://olivier.chapelle.cc/pub/hessmkl.pdf

t=0;
tic;

use_svqp = 1;  % If 0, use Matlab's quadprog (much slower)

global K; K = K_; clear K_; % Global variable to gain speed 
  
n = length(Y);   % Number of points
if size(Y,1)==1, Y=Y'; end
d = size(K,3);   % Number of kernels

if (exist('svqp')~=3) && use_svqp
  mex svqp.cpp
end;

% If C is not given, do hard margin. But if one of the kernel is
% the identity, this is equiavelent to an L2 penalization of the slacks.
if ~(exist('C'))
  C = 1;
end;

if ~(exist('max_dgap'))
  max_dgap = -1;
end;

dgap =inf;
p_obj = inf;
d_obj = inf;

% Check the kernels  
for i=1:d
  K1 = K(:,:,i);
  if abs(sum(K1(:))) > 1e-5
    %error('Kernel not centered in feature space');
  end;
  if abs(trace(K1)- n) > 1e-5
    %error('Kernel not normalized');
  end;
end;

D=ones(d,1) / d;     % Initialization: put equal weight to all kernels
alpha = zeros(n,1);

options = optimset('display','off','LargeScale','off');

% Do the optimization: 
% minimize the SVM objective function under constraints D >= 0 and sum(D)=1
% For that, compute the gradient and the Hessian of the SVM objective
% function with respect to D and make a constrained Newton step. 

iteration =0;
while 1
  iteration = iteration +1;
  %fprintf('iteration = %d\n',iteration);
  D(D<1e-10) = 0;
  [w,grad,hess,alpha] = obj_fun(D,Y,C,alpha,use_svqp); 
	
  hess = (hess+hess')/2 + 1e-5*eye(size(hess))*mean(diag(hess)); % Because of numerical errors
  fprintf('Obj fun = %f, L0 norm of the combination = %d    \n',w,sum(D>0));
  
  % Find the search direction. 
  % If the problem were unconstrained, it would be hess \ grad (Newton step).
  if use_svqp
     S = svqp(hess,-grad,-D,inf*ones(d,1),zeros(d,1),1);
     fval = 0.5*S'*hess*S +S'*grad;
  else
    [S,fval] = quadprog(hess,grad,[],[],ones(1,d),0,-D,[],[],options);
  end; 

	if max_dgap<=0
    if fval > -1e-5*w break; end; % Stop when the relative increase to the
																	% objective function is small.
	else
	  t = t + toc;				
    [dgap,p_obj,d_obj] = dualitygap(abs(alpha), Y, C, D);
	  if dgap < max_dgap, break; end;																
	  tic;
	end	
	

  % Back tracking if necessary [usually, lambda = 1]
  lambda = 1; non_convex=0;
  while 1
    w1 = obj_fun(D+lambda*S,Y,C,alpha,use_svqp);
    if w1<w break; end;
    lambda = lambda / 2;
    if lambda<1e-10 non_convex=1; break; end;
  end;
  
  %Take the step
  D = D+lambda*S;
  
  if non_convex
    warning('Convexity problem; lambda too small');
    %break;
  end;
end;
fprintf('\nNumber of MKL iterations:  %d\n', iteration); 
t=t+toc;

  
function [w,grad,hess,alpha] = obj_fun(D,Y,C,alpha0,use_svqp)
  % Compute the SVM objective function, as well as the gradient and
  % Hessian with respect to D.
  % Note that alpha depends on D.
  global K;

  n = length(Y);
  d = length(D);

  % Compute the convex combination of kernels
  K0 =  1e-8*eye(n); % Small ridge for numerical reasons
  for i=1:d
    K0 = K0 + K(:,:,i)*D(i);
  end;
  K0 = (K0+K0') / 2; % Make sure the kernel is symmetric
  
  % Standard soft margin SVM training
  if use_svqp
    alpha = svqp(K0,Y,-C*(Y==-1),C*(Y==1),alpha0,1);
  else
    options = optimset('Display','off','LargeScale','off');
    alpha = quadprog(K0,-Y,[],[],ones(1,n),0,...
                     -C*(Y==-1),C*(Y==1),alpha0,options);
  end;
    
  % For L1 SVM, we have two kind of SVs: sv1 are such that (0<alpha<=C)
  % and sv1(sv2) are the active ones (0<alpha<C)
  sv1 = find(alpha);
  sv2 = find(abs(alpha(sv1))<C);
    
  % Objective function
  w = sum(alpha.*Y) - 0.5*alpha(sv1)'*K0(sv1,sv1)*alpha(sv1);
  if nargout == 1,
    return;
  end;
  
  % Compute the gradient and the hessian of the margin
  lsv = length(sv2);
  invKsv = inv([[K0(sv1(sv2),sv1(sv2)) ones(lsv,1)]; [ones(1,lsv) 0]]);

  for i=1:d
    alK(:,i) = K(sv1,sv1,i)*alpha(sv1);
  end;
  hess = 2*alK(sv2,:)'*invKsv(1:end-1,1:end-1)*alK(sv2,:);
  grad = - alK'*alpha(sv1);
 

function [dgap,prim_obj,dual_obj] = dualitygap(alpha, Y, C, D);
  
	global K;
		
  Kw = zeros(size(K,1),size(K,2));
  for i=1:size(K,3)
	  quadterm(i) = 0.5*(alpha.*Y)'*K(:,:,i)*(alpha.*Y);
		Kw = Kw + D(i)*K(:,:,i);
	end
  ind = find((alpha>0).*(alpha<C));
	l = length(ind);
	b = 1/l*sum(Y(ind)'-(alpha.*Y)'*Kw(:,ind));
  dual_obj = sum(alpha) - max(quadterm);
  prim_obj = 0.5*(alpha.*Y)'*Kw*(alpha.*Y) + C*sum(max(0,1-(((Y.*alpha)'*Kw+b).*Y')));
  assert(dual_obj<=prim_obj);
  dgap = abs(1-dual_obj/prim_obj);	
	
	