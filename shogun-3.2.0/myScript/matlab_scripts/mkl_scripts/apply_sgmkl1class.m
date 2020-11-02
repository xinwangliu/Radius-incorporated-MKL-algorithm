function out= apply_sgmkl1class(varargin)

%apply_sgmkl1class - applies data to classifier 'classy' trained by train_mkl1class  (using C and shogun instead of nu-SVM)
%
% Synopsis:
%   s = apply_mklSVM( 'svm', svm, 'K',K)
%
% Properties:
%   svm:   trained SVM structure
%   K:  [ntr,nte, k]     testing kernel matrices
%
% Returns:
%  s:    -  real label for each test example (decision should be
%            done on the sign)


properties= propertylist2struct(varargin{:});
properties = set_defaults(properties, ...
                         'svm', [] , ...
                         'K', [], ...
                         'linX', [], ...
                         'rbfX', [], ...
                         'mkl_norm', []);
struct2workspace(properties);

if max(size(K)>0)
  svm.alphas = svm.alpha;
  if size(K,3)>1 && isfield(svm,'beta') && (~(max(size(mkl_norm))>0) || mkl_norm<inf)
    out=apply_mklSVM(svm,svm.beta,K);
  else
    out = apply_sgSVM( svm, [], 'Kts', sum(K,3) );
  end
elseif max(size(linX)>0)
  out = apply_sgSVM2('svm', svm, 'linX',linX );
else
  out = apply_sgSVM2('svm', svm, 'rbfX',rbfX );
end

out = -out;

