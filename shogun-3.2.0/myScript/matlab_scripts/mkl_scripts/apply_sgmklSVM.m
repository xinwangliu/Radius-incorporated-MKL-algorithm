function out= apply_sgmklSVM(varargin)

%apply_sgmklSVM - applies data to classifier 'classy' trained by train_mklSVM.
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
                         'X', []);
struct2workspace(properties);
svm.alphas = svm.alpha;

if max(size(K)>0)
  if size(K,3)>1 && isfield(svm,'beta') 
    out=apply_mklSVM(svm,svm.beta,K);
  else
    out = apply_sgSVM(  'K', sum(K,3),'svm', svm );
  end
else
  assert(max(max(size(X)))>0);
  %prepare shogun
  sg('clean_kernel');
  sg('clean_features', 'TRAIN');
  sg('clean_features', 'TEST');
  sg('set_kernel', 'COMBINED', 50);
  if svm.mkl_norm<inf
    sg('new_classifier', 'MKL_CLASSIFICATION');
  else
    sg('new_svm', 'LIGHT');
  end
  sg('send_command', 'use_linadd 1');
  for i=1:max(size(X,3),max(max(size(svm.kernel))))
    sg('add_features','TRAIN', svm.suppvec(:,:,min(i,size(svm.suppvec,3))));
    sg('add_features','TEST', X(:,:,min(i,size(X,3))));
    kernel_str = svm.kernel{min(i,length(svm.kernel))};
    sg('send_command', ['add_kernel 1 ' kernel_str]);
    if ~(upper(kernel_str(1))=='L')
      sg('use_linadd',0);
    end
    if length(svm.kernel_normalization{i}>0)
      sg('set_kernel_normalization', svm.kernel_normalization{i});
      if ~(upper(svm.kernel_normalization{i}(1))=='L')
        sg('use_linadd',0);
      end
    end
  end
  
  sg('set_subkernel_weights',svm.beta);
  sg('send_command', 'init_kernel TRAIN');
  sg('send_command', 'init_kernel TEST');
  sg('set_svm', svm.b, svm.alphas);
  out=sg('svm_classify');
end

