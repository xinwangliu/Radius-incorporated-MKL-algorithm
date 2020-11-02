function [K, sigma] = rbf(K, sigma, varargin)

% rbf - calculates rbf kernel matrix from linear kernel matrix
%
% Synopsis:
%   K = rbf(K, sigma)
%   [K, sigma] = rbf(K)
%
% Arguments:
%   K : linear kernel matrix
%   sigma:  rbf kernel width
%
% Properties
%   Ktrdg:     training kernel diag (optional)
%   Ktsdg:     testing kernel diag (optional)
%
% Returns:
%   rbf kernel matrix
%
% Description:
%   If no sigma is specified, sigma is chosen via a heuristic
%
% Author: Marius Kloft, Jun 2008

prop = propertylist2struct(varargin{:});

if (length(sigma)<=1)
  if ~exist('sigma') || sigma==0 || ~isfinite(sigma)
    for i=1:size(K,3)
      % heuristic taken from Guyon: Feature Selection, Chapter 20
      sigma(i) = mean(min( kern2dist(K(:,:,i)) + diag(inf*ones(1,size(K,1))) ));
      % alternative heuristic
      % tmp = triu(kern2dist(K(:,:,i)),1);
      % ssigma(i) = median(tmp(find(tmp>0)));
      prop= propertylist2struct(varargin{:});
      prop= set_defaults(prop, 'Ktrdg', diag(K(:,:,i)), 'Ktsdg', diag(K(:,:,i)) );
      [m n foo] = size(K);
      K(:,:,i) = exp( - ( prop.Ktrdg*ones(1,n) - 2*K(:,:,i) + ones(m,1)*prop.Ktsdg' ) / (2*sigma(i)^2 ));
    end
  else
    for i=1:size(K,3)
      prop= propertylist2struct(varargin{:});
      prop= set_defaults(prop, 'Ktrdg', diag(K(:,:,i)), 'Ktsdg', diag(K(:,:,i)) );
      [m n foo] = size(K);
      K(:,:,i) = exp( - ( prop.Ktrdg*ones(1,n) - 2*K(:,:,i) + ones(m,1)*prop.Ktsdg' ) / (2*sigma^2) );
    end
  end
else
  if size(K,3)==1 && length(sigma)>1
    Kaux=K;
    for i=1:length(sigma)
      [m n foo] = size(K);
      K(:,:,i) = exp( - ( 2*diag(Kaux)*ones(1,n) - 2*Kaux  ) / (2*sigma(i)^2) ),
    end
  else
    for i=1:size(K,3)
      prop= propertylist2struct(varargin{:});
      prop= set_defaults(prop, 'Ktrdg', diag(K(:,:,i)), 'Ktsdg', diag(K(:,:,i)) );
      [m n foo] = size(K);
      K(:,:,i) = exp( - ( prop.Ktrdg*ones(1,n) - 2*K(:,:,i) + ones(m,1)*prop.Ktsdg' ) / (2*sigma^2) );
    end
  end
end