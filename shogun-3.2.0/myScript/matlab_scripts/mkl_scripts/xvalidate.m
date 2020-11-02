function result = xvalidate(varargin)

% xvalidate - overloaded function for validation
%
% Synopsis:
%   result = xvalidate('K',K,'y',y,...)
%   
% Arguments:
%     none.
%
% Returns:
%    result struct
%                      
% Properties:
%    'K':                 kernel matrix/matrices
%    'X':                 data matrix (only for LibLPM)
%    'X':                 feature vectors (not yet implemented)
%    'y':                 labels.   1class: e.g. y=[ 0 0 0 ... 1 1 1 2 2 2 2 3 3 ....]
%                                        2class: e.g.     y= [-1 -1 1 1 -1 ....]
%                                                         OR y=[1 0  0 0 1 ...; 0 1 1 1 0 ...].
%    'one_class':   0 or 1.(default: 0)
%    'loss':             'auc' (default for 1-class) or '0_1' (default for 2-class)
%    'rep':             switches sampling divisions ON (and cv switched OFF), e.g. rep=100 (default).
%    'cv':               repetitions and folds. switches cross-validation ON and sampling divisions OFF,
%                          e.g. cv=[10 10]  %not yet implemented
%                          (ONLY USE EITHER rep (default) OR cv)
%    'balanced':      0 or 1. draw balanced train/val/test sets
%    n:                   training / validation / test data size,  e.g. n=[1000,500,500];  or n=1000 for 1-class.
%    roc_ub:          roc curve only evaluated in the small interval [0, roc_ub]
%                          (defaults:    0.01 for 1-class;  1 for 2-class)
%    'use_rbf':               0 or 1.  0==linear_kernel; 1==gaussian_kernel (default: 0)
%    'rbf_width':    RBF bandwidths to check, if rbf==1.  (default: sigma = mean(mindist(x_i,x_j)) )
%    'C':                 soft margin parameter (default:  C=Var(X)^(-1) ) ; 2-class only
%    'nu':               soft margin parameter (default:  nu=logspace(-2,0,3) ); 1-class only
%    'mkl_nsp':      0 or 1. turn non-sparse MKL on/off (default: 0, i.e. sparseMKL).
%    'solver':         mkl solver:   either 'internal' or 'cplex'
%    'verbosity':    -1, 0 or 1. (default: 0)
%    'mkl_norm':   p>=1
%    'elasticnet_lambda':     elasticnet MKL lambda
%    'epsilon':       optimization precison
%    'duality_gap': double-check optimization precision in terms of duality gap
%    'full_val':  sets parameters to following defaults:
%                                 sigma = logspace( 0.1*s, 10*s, 10),  s=mean(mindist(x_i,x_j)))
%                                 C = logspace( 0.1*c, 10*c, 10),  c=Var(X)^(-1)
%
%
% Description:
%   Validation can be called for one-class, two-class learning as well as multiple kernel learning.
%
% Examples:
%   1-class:
%
%   2-class:
%      X=[randn(5,50)-1  randn(5,50)+1];
%      y=[-ones(1,50) ones(1,50)];
%     result=xvalidate('K',X'*X,'y',y);
%
%   1-class MKL:
%
%   2-class MKL:
%
% train_mklSVM, train_SVM, train_kr_svdd, train_kr_qssvm, apply_mklSVM, ....
%
% Author: Marius Kloft, Jul 2008
%

X=[];  %not implemented yet

prop = propertylist2struct(varargin{:});
prop = set_defaults(prop, ...
                         'K', [], ...
                         'X', [], ...
                         'y', [], ...
                         'rep', 100, ...
                         'cv', [], ...
                         'balanced', 0, ...
                         'duality_gap', inf, ...
                         'rbf_width', [], ...
                         'mkl_norm', 1,  ...
                         'verbosity', 0, ...
                         'regression',0,...
						 'progress_bar',1,...
                         'implementation', 'LIBSVM', ...
                         'solver', 'CPLEX',...
                         'epsilon', 1e-3', ...
						 'block_norm', [], ...
						 'elasticnet_lambda', [] );
						 

prop = set_defaults(prop, 'use_rbf', size(prop.rbf_width,1)>0 );

struct2workspace(prop);
clear prop;
max_duality_gap =duality_gap;
clear duality_gap;

if exist('one_class')
  if one_class==1
    multi_class=0;
  else
    if  length(find(y==0))>0
      multi_class=1;
	else
	  multi_class=0;
	end
  end
end  
if ~exist('multi_class')&&~exist('one_class')
  a = length(find(y==0));
  b = length(find(y>1));
  if a==0
    one_class=0;
	multi_class = 0;
  elseif b>a
    one_class=0;
	multi_class = 1;
  else
    one_class=1;
	multi_class = 0;
  end
end
if multi_class
  one_class=0;
  if ~exist('mc_svm')
    mc_svm = 'SCATTER';
  end
end
if one_class
  multi_class=0;
end

if length(block_norm)>0,
  mkl_norm=block_norm;
  flag_block_norm = 1;
else
  flag_block_norm = 0;
end

if length(elasticnet_lambda)>0,
  mkl_norm=1;
  flag_elasticnet_MKL = 1;
else
  flag_elasticnet_MKL = 0;
end

if mkl_norm>1
  mkl_nsp =1;
  if mkl_norm==inf
    K=sum(K,3);
  end
else
  mkl_nsp=0;
end
K=addRidge(K);

if size(y,1)==1 && length(find(y>0))>0 &&length(find(y<0))>0
  y= abs([y+1;y-1]/2);
end

% Assertions
if ~exist('one_class')
  if  size(y,1) ==1
    if any(y==-1)
      one_class = 0;
    else
      one_class = 1;
    end
  else
    if size(y,1)==2
      one_class=0;
    else
      one_class=1;
    end
  end
end
if one_class
  if exist('nu')
    C= nu;
    isC=false;
  else
    isC = true;
  end
end
if size(K,3)<2
  mkl_norm=inf;
elseif mkl_nsp==1
  mkl_norm = 2;
else
  mkl_norm = 1;
end
if exist('ntr')
  n = ntr;
end
if ~exist('use_rbf')
  use_rbf  = ( size(rbf_width) > 0 );
end
v = kvar( mean(K,3) );
if exist('full_val')&&full_val
  if ~exist('C')&&~exist('nu')
    C = logspace( log10(0.1*v), log10(10*v), 10);
    nu = logspace(-3,0,10);
  end
else
  full_val=0;
  if ~exist('C')
    C = v;
    nu = 0.1;
  end
end
%if (length(find(abs(diff(sum(y,1)))>0.0000001))>49)
%  regression = 1;
%end
if regression
  one_class=0;
end
assert( (size(X,2)==0) || (size(K,1)==0) );
assert( ~( size(X,2)==0 & size(K,1)==0 ) );
if one_class==0&&~regression&&~multi_class
  if size(y,1) ==1
    if all(y>=0)
      y = (logical(y)-0.5)*2;
    end
    y = abs([y+1; y-1]/2);
  end
  [foo, large_class] = max( [ length(find(y(1,:)==0))  length(find(y(2,:)==0)) ] );
  if large_class==1
    y= [y(2,:);  y(1,:)];
    large_class = 2;
  end
end
if one_class==1
  if size(y,1)>1;
    y = y .* ([0:size(y,1)-1]' * ones(1,size(y,2)));
    if length(find(y==0))<length(find(y==1))
      y=not(y);
    end
    y = sum(y,1);
  elseif any(y==-1)
      y = (y+1)/2;
  end
end
if one_class
  idxn = find(y==0);
  idxa = find(y>0);
  lenATTval = 0;
  lenATTte = 0;
  for class=1:max(y(idxa))
    idx = find(y(idxa)==class);
    lenATTval = lenATTval +  length(idx(1:floor(length(idx)/2)));
    lenATTte = lenATTte + length(idx(floor(length(idx)/2)+1:end));
  end
  nAtt  = [0 lenATTval  lenATTte];
  if ~exist( 'roc_ub'), roc_ub=0.01; end
  if exist('n')
    assert(n(1)<=length(find(y==0)));
  else
    n = floor(floor(length(y)/2)*roc_ub)/roc_ub;
  end
  ntr= n(1);
  if length(n)==1
    ntot= length(idxn);
    nval = floor((floor((ntot-ntr)/2) ) *roc_ub) /roc_ub ;
    nte = floor((floor((ntot-ntr)/2) ) *roc_ub) /roc_ub ;
    %nte = floor((ntot-ntr-floor((ntot-ntr)/2) ) *roc_ub)  / roc_ub ;
    %ntr = floor((ntot-ntr-floor((ntot-ntr)/2) ) *roc_ub)  / roc_ub ;
    n= [ntr nval nte];
  end
  warning on;
  if ~(mod( n(2) , roc_ub^(-1) ) == 0)
    warning('validation set fp axis cannot be properly cutted at ub\n');
  end
  if ~(mod( n(3) , roc_ub^(-1) ) == 0)
    warning('test  set fp axis cannot be properly cutted at ub\n');
  end
  warning off;
elseif exist('n')  %2class
  if balanced&&~one_class&&~regression
    idx1 = find(y(1,:)==1);
    idx1 = idx1( randperm(length(idx1)) );
    idx2 = find(y(2,:)==1);
    idx2 = idx2( randperm(length(idx2)) );
    nn = 2*min(length(idx1),length(idx2));
  else
    nn = size(K,2);
  end
  if length(n)==1
    n(2) = floor((nn-n(1))/4) *2;
    n(3) = floor((nn-n(1))/4) *2;
  elseif length(n)==2
    n(3) = floor((nn-n(1)-n(2))/2)*2;
  else
    assert(length(n)==3);
    if balanced&&~one_class
      assert(2*min(length(idx1),length(idx2))>=sum(n));
    end
  end
else
  if balanced&&~one_class&&~regression
    idx1 = find(y(1,:)==1);
    idx1 = idx1( randperm(length(idx1)) );
    idx2 = find(y(2,:)==1);
    idx2 = idx2( randperm(length(idx2)) );
    nn = 2*min(length(idx1),length(idx2));
  else
    nn = max(size(K,2),size(X,2));
  end
  nn = floor(nn/8) *2;
  n = [2*nn, nn, nn];
end
ntr = n(1);  nval=n(2);  nte=n(3); ntotal = sum(n);
if ~exist('loss')
  if one_class
    loss = 'auc';
  else
    loss='0_1';
  end
end
if loss(1)=='A'
  loss(1)='a';
end
assert(loss(1)=='0' || loss(1)=='a');
flag_K = (length(X)==0);
flag_mkl = ( ~(size(K,1)==0) && size(K,3)>1 );

if regression
  X = empKM(K);
end
if one_class
  cv=[];
end
if any(size(cv))
  rep = cv(1);
  folds = cv(2);
  gfolds =folds;
  n = floor(size(K,2)/folds) * folds;
  balanced = 0;
else
  folds=1;
  gfolds=2;
end

%init vars
if use_rbf && (length(rbf_width)==0)
  [foo, sigma] = rbf(mean(K,3), inf);
  clear foo;
  if full_val
    rbf_width = logspace( log10( 0.1*sigma ), log10( 10*sigma ), 10 );
  else
    rbf_width = sigma;
  end
end
max_i = max(1,length(rbf_width));
max_j = length(C);
max_k = rep;
verr = -inf*ones(max_i, max_j,max_k);
terr = -inf*ones(max_i, max_j,max_k);
vauc = inf*ones(max_i, max_j,max_k);
tauc = inf*ones(max_i, max_j,max_k);

%enforce 1-d labels
if size(y,1)==2
  y = sign([y(1,:) + -y(2,:)]);
end

if progress_bar==1
  tic;
end
k=0;
c=0;
tt=0;
%REPETITIONS
if multi_class
  nr_cl = length(unique(y(:)));
  for i=1:nr_cl
    ind{i} = find(y==i-1);
    count(i) = length(ind{i});
  end
end
for kk=1:rep
  for f=1:folds
    for g=f+1:gfolds
      k = k+1;
      % 1-CLASS
	  if multi_class
	    % should be improved to BALANCED testing! Needed for many classes.
		proc_tr = n(1)/sum(n);
		proc_val = n(2)/sum(n);
		proc_te = n(3)/sum(n);
		divtr = []; divval = []; divte = [];
		for i=1:nr_cl
		  r = randperm(length(y(ind{i})));
		  divtr = [divtr ,  ind{i}(r(1:round(length(ind{i})*proc_tr)))];
		  divval = [divval ,  ind{i}(r(1+round(length(ind{i})*proc_tr):round(length(ind{i})*(proc_val+proc_tr))))];
		  divte = [divte ,  ind{i}(r(1+round(length(ind{i})*(proc_val+proc_tr)):end))];
		end
        %divtr = div(1:ntr);
        %divval = div(ntr+1 : ntr+nval);
        %divte  = div(ntr+nval+1 : ntr+nval+nte);
	  elseif one_class
        ntot=length(idxn);
        r = randperm(ntot);
        divtr = idxn(r(1:ntr));
        divval = idxn(r(ntr+[1:n(2)]));
        divte = idxn(r(ntr+n(2)+[1:n(3)]));
        for class=1:max(y(idxa))
          idx = find(y(idxa)==class);
          divval = [divval idxa(idx(1:floor(length(idx)/2)))];
          divte = [divte idxa(idx(floor(length(idx)/2)+1:end))];
        end
      % 2-CLASS
      else
        if ~exist( 'roc_ub'), roc_ub=1; end
        if ~any(size(cv))
          if balanced&&~regression
            idx1 = find(y(1,:)==1);
            idx1 = idx1( randperm(length(idx1)) );
            idx2 = find(y(2,:)==1);
            idx2 = idx2( randperm(length(idx2)) );
            l = 2*min(length(idx1),length(idx2));
            assert(l>=sum(n));
            assert( all(mod(n,2)==0) );
            divtr = [ idx1(1:ntr/2)  idx2(1:ntr/2)  ];
            divval = [ idx1(ntr/2+[1:n(2)/2])  idx2(ntr/2+[1:n(2)/2]) ];
            divte = [ idx1(ntr/2+n(2)/2+[1:n(3)/2])  idx2(ntr/2+n(2)/2+[1:n(3)/2])  ];
            divtr = divtr(randperm(length(divtr)));
            divval = divval(randperm(length(divval)));
            divte = divte(randperm(length(divte)));
          else
            div = randperm(length(y));
            divtr = div(1:ntr);
            divval = div(ntr+1 : ntr+nval);
            divte  = div(ntr+nval+1 : ntr+nval+nte);
          end
        else % cv
          if f==1 && g==2
            div = randperm(length(y));
            div = div(1:n);
          end
          len = floor(length(div)/folds);
          divte = div( (f-1)*len+1 : f*len );
          divval = div( (g-1)*len+1 : g*len );
          divtr = div( [1:(f-1)*len  f*len+1:(g-1)*len  g*len+1:folds*len] );
        end
      end
      
      %RBF: loop over bandwidth
      for i = 1: max_i
        if use_rbf
          Kaux = rbf( K, rbf_width(i) );
        elseif flag_K
          Kaux = K;
        end

        %C: loop over soft margin parameter
        for j=1:length(C)
          c = c+1;
		  if flag_K
            Ktr = Kaux(divtr,divtr,:);
            Kval = Kaux(divtr,divval,:);
            Kts = Kaux(divtr,divte,:);
            if one_class
              for b=1:size(Kaux,3)
                Kvaldg(:,:,b) = diag(Kaux(divval,divval,b));
                Ktsdg(:,:,b) = diag(Kaux(divte,divte,b));
              end
            end
		  else
		    Xtr = X(:,divtr);
			Xval = X(:,divval);
			Xte = X(:,divte);
		  end
          ytr = y(:,divtr);
          yval = y(:,divval);
          yte = y(:,divte);
          if multi_class
			if upper(mc_svm(1))=='L'
			  classy = train_mcsvm( 'K', Ktr, 'C', C(j), 'y', ytr, 'classy', 'LIBSVM_MULTICLASS');
			  outval = test_mcsvm('K', Kval, 'classy','LIBSVM_MULTICLASS', 'y', ytr);
			  outte = test_mcsvm('K', Kts, 'classy','LIBSVM_MULTICLASS', 'y', ytr);
			elseif upper(mc_svm(1))=='S' 
              classy = train_mcsvm( 'K', Ktr, 'C', C(j), 'y', ytr, 'classy', 'SCATTER');
			  outval = test_mcsvm('K', Kval, 'classy','SCATTER', 'y', ytr);
			  outte = test_mcsvm('K', Kts, 'classy','SCATTER', 'y', ytr);
			else 
			  error('invalid mc svm');
			end
            verr(i,j,k) = mean(outval(:)~=yval(:));
            terr(i,j,k) =  mean(outte(:)~=yte(:));
          elseif regression
            model = train_gfSVR(X(:,divtr), ytr, 'C', C(j));
            outval = apply_gfSVR( model, X(:,divval));
            outte = apply_gfSVR( model, X(:,divte));
            vmae(i,j,k) = mean(abs(yval-outval));
            tmae(i,j,k) = mean(abs(yte-outte));
	      % 1-CLASS
          elseif one_class
            if flag_mkl
              paraSVDD.nu=C(j);
              paraSVDD.verbosity=verbosity;
              paraSVDD.eps_mkl=epsilon;
%                if isC
%                  svm = train_sgmkl1class( 'K', Ktr, 'C', C(j), 'epsilon', epsilon, 'verbosity', verbosity, 'mkl_norm', mkl_norm );
%                  beta = svm.beta;
%                  outval = apply_sgmkl1class('svm', svm, 'K', Kval );
%                  outte = apply_sgmkl1class('svm', svm, 'K', Kts );
%                else
                [ beta, classy, resMKL ] = train_mklSVM1class( Ktr, 'nu', paraSVDD.nu, 'verbosity', verbosity,'mkl_norm',mkl_norm );
                outval = apply_mklSVDD( classy, beta, 'Kts', Kval, 'Kdg', Kvaldg );
                outte = apply_mklSVDD( classy, beta, 'Kts', Kts, 'Kdg', Ktsdg );
%              end
              [foo ,fooo, auc]=val_ROC( outval , yval, 'ub', roc_ub);
              vauc(i,j,k) = auc;
              [tp ,fp,auc]=val_ROC( outte , yte, 'ub', roc_ub);
              tauc(i,j,k) = auc;
              roc{i}{j}(k).tp = tp;
              roc{i}{j}(k).fp = fp;
              betas{i}{j}(k,:) = beta;
             else
             if isC
               classy = train_sg1classSVM( [] , ytr, 'C', C(j), 'epsilon', epsilon, 'verbosity',   verbosity, 'Ktr', Ktr );
               outval = apply_sg1classSVM(classy, [], 'Kts', Kval );
               outval = apply_sg1classSVM(classy, [], 'Kts', Kts );
             else
                if C==1
                  classy = train_kr_qssvm([], 'Ktr', Ktr );
                else
                  classy = train_kr_svdd([], 'Ktr', Ktr, 'nu',C(j) );
                end
                outval = apply_kr_svdd( classy, [], 'Kts', Kval, 'Kdg', Kvaldg );
                outte = apply_kr_svdd( classy, [], 'Kts', Kts, 'Kdg', Ktsdg );
              end
              [foo ,fooo, auc]=val_ROC( outval , yval, 'ub', roc_ub);
              vauc(i,j,k) = auc;
              [tp ,fp,auc]=val_ROC( outte , yte, 'ub', roc_ub);
              tauc(i,j,k) = auc;
              roc{i}{j}(k).tp = tp;
              roc{i}{j}(k).fp = fp;
            end
          %2-CLASS
          else
            if flag_mkl
%obsolete code
%                  if mkl_norm>1
%                    [beta, svm] = train_mklSVMnsp( Ktr , ytr, 'C', C(j), 'implementation', implementation, 'epsilon', epsilon, 'verbosity', verbosity, 'mkl_norm', mkl_norm );
%                  else
%                    [beta, svm] = train_mklSVM( Ktr , ytr, 'C', C(j), 'implementation', implementation, 'epsilon', epsilon, 'verbosity', verbosity );
%                  end
			  if flag_elasticnet_MKL==1
  			    %svm = train_lbfgsbMKL( 'K', Ktr, 'y', ytr, 'C',C(j),'verbosity',verbosity, 'mkl_norm', 1.05,'lambda',elasticnet_lambda);
				svm = train_sgmklSVM( 'K', Ktr , 'y', ytr, 'C', C(j), 'epsilon', epsilon, 'verbosity', verbosity, 'elasticnet_lambda', elasticnet_lambda,'duality_gap', max_duality_gap,'solver', solver );			  
			  elseif flag_block_norm==1
			    svm = train_lbfgsbMKL( 'K', Ktr, 'y', ytr, 'C',C(j),'verbosity',verbosity, 'mkl_norm', block_norm);
			  else
                svm = train_sgmklSVM( 'K', Ktr , 'y', ytr, 'C', C(j), 'epsilon', epsilon, 'verbosity', verbosity, 'mkl_norm', mkl_norm,'duality_gap', max_duality_gap,'solver', solver );			  
				duality_gap(i,j,k)=svm.duality_gap;
			  end
              outval = apply_sgmklSVM('svm', svm, 'K', Kval );
              outte = apply_sgmklSVM('svm', svm, 'K', Kts );
              betas{i}{j}(k,:) = svm.beta;
            elseif flag_K
              classy = train_sgmklSVM('y', ytr, 'C', C(j), 'implementation', implementation, 'epsilon', epsilon, 'verbosity',  verbosity, 'K', Ktr,'duality_gap', max_duality_gap,'mkl_norm', inf );
              outval = apply_sgmklSVM('svm', classy, 'K', Kval );
              outte = apply_sgmklSVM('svm',classy, 'K', Kts );
              duality_gap(i,j,k)=classy.duality_gap;
			else

			  classy = train_sgLPM('y', ytr, 'C', C(j),'X',Xtr);
			  outval = apply_sgLPM('X', Xval );
			  outte = apply_sgLPM('X', Xte );
            end
            if loss(1)=='0'
              verr(i,j,k) = mean(sign(outval)~=yval);
              terr(i,j,k) =  mean(sign(outte)~=yte);
            else
              if size(y,1)==1
                large_class=1;
                yval_aux=yval;
                yte_aux=yte;
                yval_aux(yval_aux==-1)=0;
                yte_aux(yte_aux==-1)=0;
              end
              [foo, fooo, auc]=val_ROC( outval , yval_aux(large_class,:) , 'ub', roc_ub);
              vauc(i,j,k) = auc;
              [tp,fp,auc]=val_ROC( outte , yte_aux(large_class,:) , 'ub', roc_ub);
              tauc(i,j,k) = auc;
              roc{i}{j}(k).tp = tp;
              roc{i}{j}(k).fp = fp;
            end
          end
		  if verbosity>-1 && progress_bar==1
            if any(cv)
              print_progress(c, max_i*max_j*max_k*folds*(folds-1)/2);
            else
              print_progress(c, max_i*max_j*max_k);%fprintf('\n');
            end
	      end
        end
      end
    end
  end
end

%Choose best model and collect results
if regression
  [foo, i] = max(max(mean(vmae,3),[],2),[],1);
  [foo, j] = max(max(mean(vmae,3),[],1),[],2);
  result.verbose.mae = mean(tmae,3);
  tmae = squeeze(tmae(i,j,:));
  result.mae= mean(tmae);
  result.std = std(tmae);
  result.stderr = std(tmae)/sqrt(rep);
elseif loss(1)=='a'  %auc
  [foo, i] = max(max(mean(vauc,3),[],2),[],1);
  [foo, j] = max(max(mean(vauc,3),[],1),[],2);
  result.verbose.auc = mean(tauc,3);
  tauc = squeeze(tauc(i,j,:));
  result.auc= mean(tauc);
  result.std = std(tauc);
  result.stderr = std(tauc)/sqrt(rep);
  result.fp = linspace(0,roc_ub, 1001);
  result.tp = mean(val_align_ROC( roc{i}{j},  result.fp ) ,1 );
elseif loss(1)=='0'
  [foo, i] = min(min(mean(verr,3),[],2),[],1);
  [foo, j] = min(min(mean(verr,3),[],1),[],2);
  result.verbose.terr = mean(terr,3);
  terr = squeeze(terr(i,j,:));
  result.err = mean(terr);
  result.stderr= std(terr)/sqrt(rep);
  result.std= std(terr);
end

if flag_mkl
  result.verbose.betas = betas;
  result.beta = mean(squeeze(betas{i}{j}),1);
end

result.rep = rep;

if one_class&&~regression&&~isC
  result.nNor = n;
  result.nAtt = nAtt;
  result.nu = C(j);
  result.nus = C;
else
  result.n = n;
  result.C = C(j);
  result.Cs = C;
end
result.epsilon = epsilon;
if use_rbf&&length(rbf_width)>0
  result.rbf_width = rbf_width(i);
end
result.rbf_widths = rbf_width;
if verbosity>-1 && progress_bar==1
  fprintf('Time: %d sec \n', round(toc));
end
result.verbose.aucx  = 'rbf';
result.verbose.aucy  = 'C';
result.date= date;
if exist('duality_gap')
  result.duality_gap=duality_gap;
end
c = clock;
result.time= [int2str_leading_zeros(c(4),2) ':' int2str_leading_zeros(c(5),2) ];
if exist('mkl_norm')
  if mkl_norm == inf
    result.beta = ones(1,size(K,3));
  end
  if flag_block_norm==1
    result.block_norm = block_norm;
  elseif flag_elasticnet_MKL==1
    result.elasticnet_lambda =  elasticnet_lambda;
  else
    result.mkl_norm = mkl_norm;
  end
end
