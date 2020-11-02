function K = kscale(K, varargin)
% scale - scales kernel matrix K (or kernel matrices K) to unit variance
%
% Synopsis:
%   K= kscale(K);
% 
% Arguments: 
%   K:     kernel matrix/matrices
%
% Return:
%   X:     kernel matrix
%
% Example:   scale(X)' *scale(X)) == kscale(X'*X)
%
% Author: Marius Kloft

warning off;

properties= propertylist2struct(varargin{:});
properties= set_defaults(properties, ...
                         'centeringANDscalekurtosis',0);

if properties.centeringANDscalekurtosis
  if ndims(K)==2
    K = kcenter(K);
    K = K./mean(diag(K).*diag(K));
  elseif ndims(K)==3
    for i=1:size(K,3)
      K(:,:,i) = kcenter(K(:,:,i));
      K(:,:,i) = K(:,:,i) ./ mean( diag(K(:,:,i)).*diag(K(:,:,i)) );
    end
  end
else
  if ndims(K)==2
    K = K ./ ( mean(diag(K))-mean(K(:)) );
  elseif ndims(K)==3
    for i=1:size(K,3)
      Kaux = K(:,:,i);
      K(:,:,i) = Kaux ./ ( mean(diag(Kaux))-mean(Kaux(:)) );
    end
  end
end