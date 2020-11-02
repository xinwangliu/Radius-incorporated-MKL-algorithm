function s = apply_sg1classSVM(C, Xtr, varargin);
%apply_sgSVM - applies data to classifier 'classy' trained by train_sgSVM.
%
% Synopsis:
%   s = apply_sg1classSVM(classy, Xts);
%   s = apply_sg1classSVM(classy, [], 'Kts', Kts);
%
% Arguments:
%   classy:	  trained SVM structure
%   Xts:	[n,d]    	real valued testing vectors
%
% Returns:
%  s:    -  real label for each test example (decision should be
%            done on the sign)
%
% Properties:
%    'Kts':   Testing kernel matrix  k(Xtr,Xts)
%
% See also: train_sgSVM
%
% Marius Kloft (2009)
% Soeren Sonnenburg (2005)
% based on work from die Guido Dornhege (2003)

properties= propertylist2struct(varargin{:});
properties= set_defaults(properties, ...
                         'Kts', [], ...
                         'Kte', [] );

if size(properties.Kts,1)==0
  properties.Kts = properties.Kte;
end
properties.Kte = 0;

sg('send_command','loglevel ERROR');
if isfield(C, 'properties') & isfield(C.properties, 'verbosity'),
  switch C.properties.verbosity,
    case 0
      %sg('set_output', '/dev/null');
	    sg('send_command','loglevel ERROR');
    case 1
      sg('send_command','loglevel WARN');
    case 2
      sg('send_command','loglevel ALL');
    otherwise
      error('Unknown value for option ''verbosity''');
  end
end

if size(C.alphas,1) > 0,
  if size(properties.Kts,1)==0
    % we need some svm object, it does not matter which one

    sg('send_command', 'new_svm LIGHT');
    sg('set_features','TRAIN', C.suppvec);
    sg('set_features', 'TEST', Xtr);
  
    sg('set_svm', C.b, C.alphas);
  
    sg('send_command', sprintf('set_kernel %s', C.kernel_string));
    sg('send_command', 'init_kernel TRAIN');
    sg('send_command', 'init_kernel TEST');

    s=sg('svm_classify');
  else
    alphas = C.alphas(:,1);
    alphas = alphas(:);
    svidx = C.alphas(:,2);
    s = C.b + alphas' * properties.Kts(svidx+1,:);
  end
else
  if size(properties.Kts,1)==0
	  s=C.b*ones(1,size(Xtr,2));
  else
    s=C.b*ones(1,size(properties.Kts,2));
  end
end

s=-s;