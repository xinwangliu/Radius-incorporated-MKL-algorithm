function v = kvar(K, nu)
% kvar  - returns variance of a kernel matrix - alternatively SVDD objectives
%
% Synopsis:
%   v = kvar(K);
%   v = kvar(K,. nu);
% 
% Arguments: 
%   K:     kernel matrix/matrices
%  nu:   SVDD nu  (default = 1)
% Return:
%   v:     vector of variances / SVDD objective values
%
% Author: Marius Kloft, Aug 2008

warning off;

if ~exist('nu')||nu==1
  for i=1:size(K,3)
    Kaux = K(:,:,i);
    v(i) = mean(diag(Kaux)) - mean(Kaux(:));
  end
else
  for i=1:size(K,3)
    for j = 1:10
      div = randperm(size(K,2));
      div=div(1:min(100,length(div)));
      C=train_kr_svdd([],'Ktr',K(div,div,i),'nu', nu);
      v(j,i) = C.objVal;
    end
  end
  v = median(v,1);
end
