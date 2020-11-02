function [D,t]=test_combination()

  % Construct a toy problem where the scale of the data is not the same
  % at different location in space. Ideally, one should combine two RBF kernels with
  % different bandwidths
  n = 50;
	C=10;
  X1 = randn(n,2);
  Y1 = sign(sum(X1.^2,2)-0.7);
  X2 = 10*randn(n,2) + repmat([10 10],n,1);;
  Y2 = sign(sum((X2-10).^2,2)-50);
  X = [X1; X2]; Y = [Y1; Y2];
  
  plot(X(Y== 1,1),X(Y== 1,2),'b+'); hold on;  
  plot(X(Y==-1,1),X(Y==-1,2),'ro'); hold off; pause(1); drawnow;

  % The base kernels are RBF kernels with various bandwidths.
  sig = 2.^[-2:4];
  for i=1:length(sig)
    K(:,:,i) = compute_kernel(X,X,sig(i));
  end;
  % And one of them is a ridge (to penalize the training errors).
  K(:,:,end+1) = eye(2*n); 
  
  % Normalize the kernels (centering in feature space and variance 1)
  %K = normalize_kernel(K);
  
	%Compute the linear combination
  [D,t,iteration,p_obj,d_obj,dgap] = kernel_combination(K,Y,C);

  % Plot the result
  figure;
  bar([1:length(sig) length(sig)+2],D);
  lab = num2str(sig');
  lab(end+1,1:5)='ridge';
  set(gca,'XTickLabel',lab);
  xlabel('Bandwidth');
  ylabel('Weight');
  
  fprintf(['The convex combination puts more weight on small and large ' ...
           'bandwidths, according to the structure of the problem.\n']);
  
function K = normalize_kernel(K)
 % Normalize the kernels (centering in feature space and variance 1)
   n = size(K,1);
  for i=1:size(K,3)
    K1 = K(:,:,i);
    S = mean(K1,2);
    K1 = K1 -  S*ones(1,n) - ones(n,1)*S' + mean(S);
    K1 = K1 / mean(diag(K1));
    K(:,:,i) = K1;
  end;
    
  
function K = compute_kernel(X,Y,sig)
  % Compute the RBF kernel matrix
  if length(sig)==1
    X = X / sig;
    Y = Y / sig;
  else
    X = X ./ repmat(sig',size(X,1),1);
    Y = Y ./ repmat(sig',size(Y,1),1);
  end;

  normX = sum(X.^2,2);
  normY = sum(Y.^2,2);

  K = exp(-0.5*(repmat(normX ,1,size(Y,1)) + ...
                repmat(normY',size(X,1),1) - 2*X*Y'));
  
