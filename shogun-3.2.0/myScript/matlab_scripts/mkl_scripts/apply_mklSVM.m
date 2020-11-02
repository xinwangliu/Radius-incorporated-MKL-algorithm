function [s,K] = apply_mklSVM(classy, beta , Kts)

%apply_mklSVM - applies data to classifier 'classy' trained by train_mklSVM.
%
% Synopsis:
%   s = apply_mklSVM(svm, beta , Kts)
%
% Arguments:
%   svm:   trained SVM structure
%   beta:     mkl kernel weights
%   Kts:  [n,n, k]     testing kernel matrices
%
% Returns:
%  s:    -  real label for each test example (decision should be
%            done on the sign)

[d, n, k] = size(Kts);
beta=beta(:);

% Tensor Multiplikation entlang der dritten Stufe
K = reshape( reshape( Kts, d*n, k) * beta  ,  d,n,1);

if ~isfield(classy,'alphas')
  classy.alphas = classy.alpha;
end
s = apply_sgSVM('svm',classy,'K', K );