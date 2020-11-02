
function [ res, info ] = trainAndTest( param, y, Ks, perms, info0 );



% === Init

N = length( y );
idxTst = param.idxTst;
idxTrn = param.idxTrn;
nofTrn = length( idxTrn );
nofKernels = size( Ks, 3 );
classes = unique( y );
nofClasses = length( classes );
ASSERT( param.foldVal == 0 );

if( param.mklNorm < inf )
  sg( 'new_classifier', 'MKL_CLASSIFICATION' );
  beta0 = [];
  sg( 'set_subkernel_weights', ones(nofKernels,1)' );
else
  fprintf( 'switching off MKL; ' );
  sg( 'new_classifier', 'SVMLIGHT' );
  if( isfield(param,'fixBetas') )
    beta0 = param.fixBetas;
    fprintf( 'using specified betas\n' );
  else
    beta0 = ones( nofKernels, nofClasses );
    fprintf( 'using uniform betas\n' );
  end;
end;



% === Train

ASSERT( strcmp( param.mcMode, '1vsRest' ) );
alphas = repmat( nan, nofTrn, nofClasses );
biases = repmat( nan, 1,      nofClasses );
betas = repmat( nan, nofKernels, nofClasses );
objs = repmat( nan, 3, nofClasses );

for( m = 1:nofClasses )

  % --- kernel weights
  if( ~isempty(beta0) )
    sg( 'set_subkernel_weights', beta0(:,m)' );
  end;

  % --- 1 vs rest
  y_m = 2*( y == classes(m) ) - 1;
  sg( 'set_labels', 'TRAIN', y_m(idxTrn)' );
  t0 = cputime;
  sg( 'train_classifier' );
  trainTimes(m) = cputime - t0;
  
  % --- retreive SVM
  [ bias, a2 ] = sg( 'get_svm' );
  alpha = zeros( nofTrn, 1 );
  alpha( a2(:,2) + 1 ) = a2(:,1);
  beta = sg( 'get_subkernel_weights' );
  
  % --- store results
  biases( :, m ) = bias;
  alphas( :, m ) = alpha;
  betas ( :, m ) = beta;
  if(  param.mklNorm < inf )
    objs( 1, m ) = -sg( 'compute_mkl_primal_objective' );
    objs( 2, m ) = -sg( 'compute_mkl_dual_objective' );
  else
    objs( 1, m ) = -sg( 'compute_svm_primal_objective' );
    objs( 2, m ) = -sg( 'compute_svm_dual_objective' );
  end;
  objs( 3, m ) = -sg( 'compute_svm_dual_objective' );
  
end;
fprintf( 'done\n' );

gaps = objs(1,:) - objs(2,:);



% === Test

% --- myself
outs = zeros( N, nofClasses );
for( k = 1:nofKernels )
  if( all( betas(k,:) == 0 ) )
    continue;
  end;
  outs = outs + Ks(:,idxTrn,k) * sparse(alphas) .* repmat(betas(k,:),N,1);
end;
outs = outs + repmat( biases, N, 1 );

% --- DEBUG: shogun
if( 0 )
  % - data
  sg( 'set_kernel', 'COMBINED', param.cacheSize );
  for( k = 1:nofKernels )
    K = Ks( idxTrn, :, k );
    sg( 'add_kernel', 1, 'CUSTOM', K, 'FULL' );
  end;
  % -- predict
  for( m = 1:nofClasses )
    y_m = 2*( y == classes(m) ) - 1;
    sg( 'set_labels', 'TEST', y_m' );
    [ is, js, as ] = find( alphas( :, m ) );
    sg( 'set_svm', biases(1,m), [ as, is-1 ] );
    sg( 'set_subkernel_weights', betas(:,m)' );
    out = sg( 'classify' )';
    max( abs( outs(:,m) - out ) )
  end;
end;



% === Return Results

[ dummy, pred ] = max( outs, [], 2 );
pred = classes( pred );
[ confmat, accuracy, tp, fp, fn, tn, precision, recall, F, MCC ] = getperf( y(idxTst), pred(idxTst), classes );

res = [];
res.F1 = F;
res.MCC = MCC;
res.acc = accuracy;
res.y = y( idxTst );
res.pred = pred;

info = info0;
info.outs = outs;
info.biases = biases;
info.alphas = alphas;
info.betas = betas;
info.objs = objs;
info.confMat = confmat;
info.trainTimes = trainTimes;

