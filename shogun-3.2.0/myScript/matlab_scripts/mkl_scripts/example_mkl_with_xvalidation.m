function result = mkl_exp_example

% Parameters
rep = 5;
n = [10 10 10];    % ntrain / nval / ntest
C = logspace(-1,1,3);	 % this is the SVM soft margin parameter;  has to be tuned per cross validation
						 % but here we simply set it to C=100, enforcing hard margin;
mkl_norm = [1, 1.333, 2, 4, inf]; 
						 % this is the crucial sparsity parameter
                         % can be any real number with  1<=mkl_norm<=infty
                         % e.g.  1-norm gives sparse solutions
                         % e.g. infty-norm:  is the regular SVM on the unweighted sum kernel
                         % e.g. 2-norm is a non-sparse MKL
loss = '0_1';
	

	
% generate some Toy data
X = [1+randn(2,20), -1+randn(2,20)];
y = [ones(1,20) -ones(1,20)];

% plot the data
figure(1); clf; hold on;
plot( X(1,1:20), X(2,1:20),  'rx');
plot( X(1,21:40), X(2,21:40), 'bo');

% employ some kernels;  you can use any kernel matrices /views here!
kernel_matrix1 = rbf( X'*X, 1);   % here we use a rbf kernel
kernel_matrix2 = X'*X;;              % and  a linear kernel
K(:,:,1)= kernel_matrix1;
K(:,:,2)= kernel_matrix2;


result = mkl_xvalidate('K', K,'y',y, 'rep', rep, 'n', n, 'C', C, 'mkl_norm', mkl_norm, 'loss', loss);


