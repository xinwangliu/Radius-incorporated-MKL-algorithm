clear;
clc
path= '/home/staff/x/xliu/Desktop/myMCMKL/KernelCode/shogun-3.2.0/';
addpath(genpath(path));
% example script for Nico

% generate some Toy data
X = [1+randn(2,20), -1+randn(2,20)];
y = [ones(1,20) -ones(1,20)];

% plot the data
figure(1); clf; hold on;
plot( X(1,1:20), X(2,1:20),  'rx');
plot( X(1,21:40), X(2,21:40), 'bo');

%train/test divisions
divtr = [1:10, 21:30];
divte = [11:20, 31:40];

% employ some kernels;  you can use any kernel matrices /views here!
kernel_matrix1 = rbf( X'*X, 1);   % here we use a rbf kernel
kernel_matrix2 = X'*X;             % and  a linear kernel
K(:,:,1)= kernel_matrix1;
K(:,:,2)= kernel_matrix2;

% MKL  Parameters
C = 100; % this is the SVM soft margin parameter;  has to be tuned per cross validation
               % but here we simply set it to C=100, enforcing hard margin;
mkl_norm = 2;  %this is the crucial sparsity parameter
                         % can be any real number with  1<=mkl_norm<=infty
                         % e.g.  1-norm gives sparse solutions
                         % e.g. infty-norm:  is the regular SVM on the unweighted sum kernel
                         % e.g. 2-norm is a non-sparse MKL


% train an MKL SVM with shogun
svm = train_sgmklSVM( 'K', K(divtr,divtr,:), 'y', y(divtr), 'mkl_norm', mkl_norm);

% calculate the test outputs, that is the distance of a test point from the separating hyperplane
% in this toy example here we just use the training data for testing
test_outputs = apply_sgmklSVM( 'K', K(divtr,divte,:), 'svm', svm);

test_error  = mean( y(divte) ~= sign(test_outputs) );
fprintf('test_error = %d%%\n', 100*test_error);


