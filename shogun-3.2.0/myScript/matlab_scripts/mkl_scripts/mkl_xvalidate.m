function result = mkl_xvalidate( varargin )

% mkl_xvalidate - validates p norm MKL and elastic net
%
% Synopsis:
%    res = train_mklSVM( properties)
%    e.g.,  res = mkl_validate('K', K,'y',y, 'rep', rep, 'n', n);
%
% Returns:
%    res:  results of validation
%
% Properties:
%        X:    data   ||   K:   kernel matrix  n x n x k
%        y:    labels  1 x n    or  2 xn
%        n:   e.g.  n=[100 200 400]   (optional)
%        loss:  'auc'  ||  '0_1'
%        mkl_norm:  1<mkl_norm<=inf   (optional)
%        elasticnet_lambda:   elastic net parameter,   0<elasticnet_lambda<1  (optional)
%        rbf_width:      rbf kernel widths  (e.g. width=logspace(-1,1,5))  (optional
%


properties= propertylist2struct(varargin{:});
properties = set_defaults(properties, ...
						'X', [],...
                        'y', [], ...
                        'K', [], ...
						'C', logspace(-2,2,5),...
						'loss', '0_1', ...
						'rep', 100, ...
                        'mkl_norm', [], ...
                        'rbf_width', [],...
					    'elasticnet_lambda', []);
struct2workspace(properties);

if length(X)>0
  assert(length(rbf_width)>0);
end
if length(mkl_norm)==0
  assert(length(elasticnet_lambda)>0);
else
  ind = find(mkl_norm<1);
  notind = setdiff([1:length(mkl_norm)], ind);
  if length(ind)>0
    elasticnet_lambda = mkl_norm(ind);
	mkl_norm = mkl_norm(notind);
  end
end

for i=1:length(mkl_norm)+length(elasticnet_lambda)
  if length(X)>0
    if i<=length(mkl_norm)
	  if exist('n')
	    res = xvalidate('X',X, 'y',y, 'rbf_width',rbf_width, 'n',n, 'C',C, 'mkl_norm',mkl_norm(i), 'loss', loss, 'rep', rep);
	  else 
	    res = xvalidate('X',X, 'y',y, 'rbf_width',rbf_width, 'C',C, 'mkl_norm',mkl_norm(i), 'loss', loss, 'rep', rep);
      end
	else
	  if exist('n')
	    res = xvalidate('X',X, 'y',y, 'rbf_width',rbf_width, 'n',n, 'C',C, 'elasticnet_lambda',elasticnet_lambda(i-length(mkl_norm)) , 'loss', loss, 'rep',rep);
	  else
	    res = xvalidate('X',X, 'y',y, 'rbf_width',rbf_width, 'C',C, 'elasticnet_lambda',elasticnet_lambda(i-length(mkl_norm)) , 'loss', loss, 'rep',rep);
	  end
	end
  else  
    if i<=length(mkl_norm)
	  if exist('n')
	    res = xvalidate('K',K, 'y',y, 'n',n, 'C',C, 'mkl_norm',mkl_norm(i), 'loss', loss, 'rep', rep);
	  else
	    res = xvalidate('K',K, 'y',y,  'C',C, 'mkl_norm',mkl_norm(i), 'loss', loss, 'rep', rep);
	  end
	else
	  if exist('n')
	    res = xvalidate('K',K, 'y',y, 'n',n, 'C',C, 'elasticnet_lambda',elasticnet_lambda(i-length(mkl_norm)), 'loss', loss, 'rep', rep);
	  else
	    res = xvalidate('K',K, 'y',y, 'C',C, 'elasticnet_lambda',elasticnet_lambda(i-length(mkl_norm)), 'loss', loss, 'rep', rep);
	  end
	end
  end
  if upper(loss(1))=='A'
    result.auc(i) = res.auc;
  else
    result.err(i) = res.err;  
  end	
  result.std(i) = res.std;
  result.stderr(i) = res.stderr;
  result.beta(i,:) = res.beta;
  result.C(i) = res.C;
end

result.mkl_norm = [mkl_norm(:)'  elasticnet_lambda(:)'];
result.rep = rep;
result.n = res.n;
result.Cs = C;
  
  