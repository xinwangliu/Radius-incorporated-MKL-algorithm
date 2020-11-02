function gap=condor_bioC_check_duality_gap(divtr, alphas, b, betas,mkl_norm)

  %load data according to training splits
  warning off;
  la = load('hg16');
  lb = load('hg16_val');
  X = [la.X, lb.X];
  y = [la.y, lb.y];
  y = [y(1,:)-0.5]*2;
  Xtr = X(:,divtr);
  ytr = y(:,divtr);

  sg('clean_kernel');
  sg('clean_features', 'TRAIN');
  sg('clean_features', 'TEST');
  sg('set_kernel', 'COMBINED', 50);
  sg('set_labels','TRAIN', ytr );
  sg('send_command', 'use_linadd 1');
  sg('use_mkl', 1);
  sg('mkl_parameters', 10^(-3), 0, mkl_norm);

%%%%%%%%%  add kernels

  %----------------WD_KERNEL--------------------------------------------
  X =Xtr;
  y = ytr;
  
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
  
  %sg('use_batch_computation', 1);
  
  sg('add_features', 'TRAIN', X,'DNA');
  sg('add_kernel',  1,'WEIGHTEDDEGREEPOS2', 'CHAR', 10, order, max_mismatch, len, shifts);
  %sg('add_features', 'TEST', X,'DNA');


  %----------PROMOTER--------------------------------------------------

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
 % sg('add_features', 'TEST', X, 'DNA');
 % sg('convert', 'TEST', 'STRING', 'CHAR', 'STRING', 'WORD', order, order-1);


  %------------EXON------------------------------------------------

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
  %sg('add_features', 'TEST', X, 'DNA');
  %sg('convert', 'TEST', 'STRING', 'CHAR', 'STRING', 'WORD', order, order-1);

    
  
  
  %---------------ANGLE---------------------------------------------
  
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
  %sg('add_features','TEST', X_ang);
  sg('add_kernel', 1,'LINEAR', 'REAL', cache, scale);
  %sg('add_kernel', 1,'POLY', 'REAL', cache, 1);
  
  
  %----------------ENERGY--------------------------------------------
  
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
%  sg('add_features','TEST', X_ene);
  sg('add_kernel', 1,'LINEAR', 'REAL', cache, scale);

%%%%%%%%%%% Ende compute Kernels

  % init Ktrain and compute duality gap
  sg('send_command', 'init_kernel TRAIN');
  %sg('send_command', 'init_kernel TEST');
  sg('send_command', 'new_svm LIGHT');
  sg('set_svm', b, alphas);
  sg('set_subkernel_weights',betas);
  %out=sg('svm_classify');
  keyboard;
  lb=sg('compute_mkl_dual_objective');
  ub=sg('compute_mkl_primal_objective');
  gap = abs(1-lb/ub);
  keyboard;

