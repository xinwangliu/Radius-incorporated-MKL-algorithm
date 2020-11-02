function K = addRidge(K, epsilon)
% addRigde - adds small ridge to diagonal of kernel matrix
%
% Synopsis:
%   K= addRidge(K);
%   K= addRidge(K, epsilon);
% 
% Arguments: 
%   K:             kernel matrix
%   epsilon:  ridge constant  (default = 10^(-10))
%                         the ridge is calculated as ridge constant *  mean(diag(K))
%
% Return:
%   K:            kernel matrix
%
%
% Author: Marius Kloft

warning off;

if ~exist('epsilon')
  epsRidge = 10^(-10);
else
  epsRidge = epsilon;
end

if ndims(K)==2
  K = K + epsRidge * mean(diag(K)) * speye(size(K,2));
elseif ndims(K)==3
  for k=1:size(K,3)
    K(:,:,k) = K(:,:,k) + epsRidge * mean(diag(K(:,:,k))) * speye(size(K,2));
  end
end
