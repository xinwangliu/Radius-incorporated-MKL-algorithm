%vars:   rep ntr divtr

if mkl_norm==1
  sg('set_solver','CPLEX');
else
  sg('set_solver','CPLEX');
end

cd('~/condor');

dirname = 'condor_bioC';
mkdir(dirname);
cd(dirname);

n= ntr;
filename = ['nr_test=', int2str(nr_test),'ntr=', int2str(ntr), 'mkl_norm=', sprintf('%0.4g',mkl_norm), ...
                      'C=',  sprintf('%0.4g',C), 'rep=', int2str_leading_zeros(rep,2), '.mat'];

if exist(filename, 'file')
  printf('Skipping saving\n');
else
  n_max = 93550;

  nr_kern = [1 2 3 4 5];
  result.nr_kern=nr_kern;
%  sg('loglevel', 'ALL');
  sg('echo', 'ON');
  %sg('progress', 'ON');
  sg('send_command','loglevel ALL');
  implementation = 'LIBSVM';
  %C = 2;
  verbosity = 0;
  mkl_eps = 10^(-3);
  svm_eps = 10^(-3);
  
  la = load('hg16');
  lb = load('hg16_val');
  X = [la.X, lb.X];
  y = [la.y, lb.y];
  
  y = [y(1,:)-0.5]*2;
   
  Xtr = X(:,divtr);
  ytr = y(:,divtr);
  
  clear X y;
  
  tic;
  
  cache=10;
  
  
  %------------------------------------------------------------
  %----------TRAINING --------------------------------------------------
  %------------------------------------------------------------
  
  %Init Kernels for Training
  
  sg('clean_kernel');
  sg('clean_features', 'TRAIN');
  sg('clean_features', 'TEST');
  sg('set_labels','TRAIN', ytr );         % set the labels
  sg('use_linadd', 1);
  sg('svm_epsilon', svm_eps);
  if mkl_norm>0
    sg('use_mkl', 1);
    sg('mkl_parameters', mkl_eps, 0, mkl_norm);
    sg('set_kernel', 'COMBINED', 0);
  else
    sg('use_mkl', 0);
    sg('set_kernel', 'COMBINED', 1);
  end
  
  %------------------------------------------------------------
  
  %----------------WD_KERNEL--------------------------------------------
  X =Xtr;
  y = ytr;
  
  if any(nr_kern ==1)
  %Spectrum kernel parameters
  order=24;
  shift=32 ;
  max_mismatch=0;
  cache=10;
  normalization='FULL'; %NO,SQRT,LEN,SQLEN,FULL
  
  X = X(1200-70:1200+70,:);
  y = ones(1,size(X,2));
  len =size(X,1);
  x=shift*rand(1,len);
  shifts = int32(floor(x(end:-1:1)));
  
  %train svm
  sg('use_batch_computation', 1);
  
  sg('add_features', 'TRAIN', X,'DNA');
  sg('add_kernel',  1,'WEIGHTEDDEGREEPOS2', 'CHAR', 10, order, max_mismatch, len, shifts);
  end
  
  %----------PROMOTER--------------------------------------------------
  if any(nr_kern ==2)
  
  X = Xtr;
  y = ytr;
  
  %Spectrum kernel parameters
  order=4;
  cache=10;
  use_sign=1;
  normalization='FULL'; %NO,SQRT,LEN,SQLEN,FULL
  
  X = X(1200-600:1200+0,:);
  y = ones(1,size(X,2));
  
  %train svm
  
  sg('add_features', 'TRAIN', X, 'DNA');
  sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'WORD', order, order-1);
  sg('add_preproc', 'SORTWORDSTRING');
  sg('attach_preproc', 'TRAIN');
  sg('add_kernel', 1,  'COMMSTRING', 'WORD', cache, use_sign, normalization);
  end
  
  %------------EXON------------------------------------------------
  if any(nr_kern ==3)
  X=Xtr;
  y = ytr;
  
  %Spectrum kernel parameters
  order=4;
  cache=10;
  use_sign=1;
  normalization='FULL'; %NO,SQRT,LEN,SQLEN,FULL
  
  X = X(1200-0:1200+900,:);
  y = ones(1,size(X,2));
  
  %train svm
  
  sg('add_features', 'TRAIN', X, 'DNA');
  sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'WORD', order, order-1);
  sg('add_preproc', 'SORTWORDSTRING');
  sg('attach_preproc', 'TRAIN');
  sg('add_kernel', 1, 'COMMSTRING', 'WORD', cache, use_sign, normalization);
  end
  
  %---------------ANGLE---------------------------------------------
  if any( nr_kern ==4)
  X = Xtr;
  y = ytr;
  
  select4 = [600:1100];
  winLen4 = 70;
  
  meanAngle = 34.1;
  par.step4 = 20;
  
  assert( length(select4) - par.step4 >= winLen4 );
  X_ang = getSmoothedDinuc( 'angle',  X(select4,:), winLen4, par.step4 ) - meanAngle;
  assert( sum( ~ isfinite( X_ang(:) ) ) == 0 );
  norms = sqrt( sum( X_ang.^2, 1 ) );
  X_ang = X_ang .* repmat( 1./norms, size(X_ang,1), 1 );
  
  scale = 1;
  sg('add_features','TRAIN', X_ang);
  sg('add_kernel', 1,'LINEAR', 'REAL', cache, scale);
  %sg('add_kernel', 1,'POLY', 'REAL', cache, 1);
  end
  
  %----------------ENERGY--------------------------------------------
  if any(nr_kern ==5)
  X = Xtr;
  y = ytr;
  
  select5 = [600:1100];
  winLen5 = 70;
  
  meanAngle = 34.1;
  meanEnergy = -7.77;
  par.step5 = 20;
  
  assert( length(select5) - par.step5 >= winLen5 );
  X_ene = getSmoothedDinuc( 'energy', X(select5,:), winLen5, par.step5 ) - meanEnergy;
  assert( sum( ~ isfinite( X_ene(:) ) ) == 0 ); 
  norms = sqrt( sum( X_ene.^2, 1 ) );
  X_ene = X_ene .* repmat( 1./norms, size(X_ene,1), 1 );
  
  scale = 1;
  sg('add_features','TRAIN', X_ene);
  sg('add_kernel', 1,'LINEAR', 'REAL', cache, scale);
  %sg('add_kernel', 1,'POLY', 'REAL', cache, 1);
  end
  
  %------------------------------------------------------------
  
  % MKL Training
  sg('c', C);
  sg('new_svm', 'LIGHT');
  sg('init_kernel', 'TRAIN');
  %sg('init_kernel_optimization');
  sg('train_classifier');
  [b,alphas]=sg('get_svm') ;
  if length(nr_kern>1)
    result.w = sg('get_subkernel_weights');
  end
  
  if ntr<100
    result.Ktr = sg('get_kernel_matrix');
  end
  
  %------------------------------------------------------------
  %-----VALIDATION-------------------------------------------------------
  %------------------------------------------------------------
  
  la = load('hg16');
  lb = load('hg16_val');
  X = [la.X, lb.X];
  y = [la.y, lb.y];
  y = [y(1,:)-0.5]*2;

  Xte = X(:,divval);
  yte = y(:,divval);
  
  clear X y;
  
  %------------------------------------------------------------
  
  fprintf('*****************VAL VAL VAL*****************');
  
  sg('clean_features', 'TEST');
  
  if any(nr_kern ==1)
    order=24;
    shift=32 ;
    max_mismatch=0;
    cache=10;
    normalization='FULL'; %NO,SQRT,LEN,SQLEN,FULL
    X = Xte(1200-70:1200+70,:);
    y = ones(1,size(X,2));
    len =size(X,1);
    x=shift*rand(1,len);
    shifts = int32(floor(x(end:-1:1)));
    sg('add_features','TEST', X, 'DNA');
  end
  
  if any(nr_kern ==2)
  order=4;
  cache=10;
  use_sign=1;
  normalization='FULL'; %NO,SQRT,LEN,SQLEN,FULL
  X = Xte(1200-600:1200+0,:);
  y = ones(1,size(X,2));
  
  sg('add_features', 'TEST', X, 'DNA');
  sg('convert', 'TEST', 'STRING', 'CHAR', 'STRING', 'WORD', order, order-1);
  %sg('add_preproc', 'SORTWORDSTRING');
  end
  
  if any(nr_kern ==3)
  order=4;
  cache=10;
  use_sign=1;
  normalization='FULL'; %NO,SQRT,LEN,SQLEN,FULL
  X = Xte(1200-0:1200+900,:);
  y = ones(1,size(X,2));
  
  sg('add_features', 'TEST', X, 'DNA');
  sg('convert', 'TEST', 'STRING', 'CHAR', 'STRING', 'WORD', order, order-1);
  %sg('add_preproc', 'SORTWORDSTRING');
  end
  
  if any(nr_kern ==4)
  X = Xte;
  select4 = [600:1100];
  winLen4 = 70;
  meanAngle = 34.1;
  par.step4 = 20;
  assert( length(select4) - par.step4 >= winLen4 );
  X = getSmoothedDinuc( 'angle',  X(select4,:), winLen4, par.step4 ) - meanAngle;
  assert( sum( ~ isfinite( X(:) ) ) == 0 );
  norms = sqrt( sum( X.^2, 1 ) );
  X = X .* repmat( 1./norms, size(X,1), 1 );
  sg('add_features','TEST', X);
  end
  
  if any(nr_kern ==5)
  X = Xte;
  select5 = [600:1100];
  winLen5 = 70;
  meanAngle = 34.1;
  meanEnergy = -7.77;
  par.step5 = 20;
  assert( length(select5) - par.step5 >= winLen5 );
  X = getSmoothedDinuc( 'energy', X(select5,:), winLen5, par.step5 ) - meanEnergy;
  assert( sum( ~ isfinite( X(:) ) ) == 0 );
  norms = sqrt( sum( X.^2, 1 ) );
  X = X .* repmat( 1./norms, size(X,1), 1 );
  sg('add_features','TEST', X);
  end
  
  %sg('add_preproc', 'SORTWORDSTRING');
  sg('attach_preproc', 'TEST');
  yte = [yte(1,:)-0.5]*2;
  sg('set_labels','TEST', yte);
  sg('init_kernel', 'TEST');
  out=sg('classify');
  yte = 0.5*(yte+1);
  result.auPRC_val = auprc(out, yte);
  [tp,fp,auc] = val_ROC(out, yte);
  result.auROC_val  = auc;
  
  %------------------------------------------------------------
  
  result.outval = out;
  result.yval = yte;


%---------TEST TEST
  la = load('hg16');
  lb = load('hg16_val');

  X = [la.X, lb.X];
  y = [la.y, lb.y];
  y = [y(1,:)-0.5]*2;
  
  Xte = X(:,divte);
  yte = y(:,divte);
  
  %------------------------------------------------------------
  
  fprintf('*****************TEST TEST TEST*****************');
  
  sg('clean_features', 'TEST');
  
  if any(nr_kern ==1)
    order=24;
    shift=32 ;
    max_mismatch=0;
    cache=100;
    normalization='FULL'; %NO,SQRT,LEN,SQLEN,FULL
    X = Xte(1200-70:1200+70,:);
    y = ones(1,size(X,2));
    len =size(X,1);
    x=shift*rand(1,len);
    shifts = int32(floor(x(end:-1:1)));
    sg('add_features','TEST', X, 'DNA');
  end
  
  if any(nr_kern ==2)
  order=4;
  cache=10;
  use_sign=1;
  normalization='FULL'; %NO,SQRT,LEN,SQLEN,FULL
  X = Xte(1200-600:1200+0,:);
  y = ones(1,size(X,2));
  
  sg('add_features', 'TEST', X, 'DNA');
  sg('convert', 'TEST', 'STRING', 'CHAR', 'STRING', 'WORD', order, order-1);
  %sg('add_preproc', 'SORTWORDSTRING');
  end
  
  if any(nr_kern ==3)
  order=4;
  cache=10;
  use_sign=1;
  normalization='FULL'; %NO,SQRT,LEN,SQLEN,FULL
  X = Xte(1200-0:1200+900,:);
  y = ones(1,size(X,2));
  
  sg('add_features', 'TEST', X, 'DNA');
  sg('convert', 'TEST', 'STRING', 'CHAR', 'STRING', 'WORD', order, order-1);
  %sg('add_preproc', 'SORTWORDSTRING');
  end
  
  if any(nr_kern ==4)
  X = Xte;
  select4 = [600:1100];
  winLen4 = 70;
  meanAngle = 34.1;
  par.step4 = 20;
  assert( length(select4) - par.step4 >= winLen4 );
  X = getSmoothedDinuc( 'angle',  X(select4,:), winLen4, par.step4 ) - meanAngle;
  assert( sum( ~ isfinite( X(:) ) ) == 0 );
  norms = sqrt( sum( X.^2, 1 ) );
  X = X .* repmat( 1./norms, size(X,1), 1 );
  sg('add_features','TEST', X);
  end
  
  if any(nr_kern ==5)
  X = Xte;
  select5 = [600:1100];
  winLen5 = 70;
  meanAngle = 34.1;
  meanEnergy = -7.77;
  par.step5 = 20;
  assert( length(select5) - par.step5 >= winLen5 );
  X = getSmoothedDinuc( 'energy', X(select5,:), winLen5, par.step5 ) - meanEnergy;
  assert( sum( ~ isfinite( X(:) ) ) == 0 );
  norms = sqrt( sum( X.^2, 1 ) );
  X = X .* repmat( 1./norms, size(X,1), 1 );
  sg('add_features','TEST', X);
  end
  
  %sg('add_preproc', 'SORTWORDSTRING');
  sg('attach_preproc', 'TEST');
  yte = [yte(1,:)-0.5]*2;
  sg('set_labels','TEST', yte);
  sg('init_kernel', 'TEST');
  out=sg('classify');
  yte = 0.5*(yte+1);
  result.auPRC_te = auprc(out, yte);
  [tp,fp,auc] = val_ROC(out, yte);
  result.auROC_te  = auc;
  
  %------------------------------------------------------------
  
  fprintf('Time: %d sec \n', round(toc));
  result.time = round(toc);
  result.outte = out;
  result.yte = yte;
  result.eps = mkl_eps;

  
  epsi = result.eps;
  struct2workspace(result);
  save(filename, 'auPRC_val', 'auROC_val', 'auPRC_te', 'auROC_te', 'outval','yval','outte', 'yte','time','epsi','w', 'divval', 'divtr', 'divte','b','alphas');

end

