
function [ info0 ] = prepForMkl( param, y, Ks, perms );



% === Split Data

N = length( y );
classes = unique( y );
nofClasses = length( classes );
numPerm = param.numPerm;

% --- test set (tst)
nofTst   = round( N * param.fracTst );
nofNoTst = N - nofTst;
idxTst   = perms( numPerm, 1:nofTst );
idxNoTst = perms( numPerm, (nofTst+1):end );
ASSERT( length(unique(y(idxNoTst))) == nofClasses );

% --- validation sets (val)
nofFoldsVal = param.nofFoldsVal;
foldVal = param.foldVal;
if( foldVal > 0 )
  % -- CV
  inFoldVal = mod( (1:nofNoTst)-1, nofFoldsVal ) + 1;
  idxVal = idxNoTst( inFoldVal == foldVal );
  idxTrn = idxNoTst( inFoldVal ~= foldVal );
  ASSERT( length(unique(y(idxTrn))) == nofClasses );
else
  % -- no CV
  ASSERT( foldVal == 0 );
  idxVal = idxTst;
  idxTrn = idxNoTst;
end;

% --- set sizes
nofVal = length( idxVal );
nofTrn = length( idxTrn );



% === Prepare Shogun

sg( 'clean_kernel' );
sg( 'clean_features', 'TRAIN' );

% --- init SVM
sg( 'svm_use_bias', 1 );
%sg( 'svm_use_bias', 0 );
sg( 'svm_epsilon', param.svmEps );
% sg( 'c', param.svmC );
%sg( 'loglevel', 'ALL' );
%sg( 'loglevel', 'INFO' );

% --- init MKL
% sg( 'use_mkl', 1 );
% sg( 'set_solver', param.mklSolver );
% sg( 'mkl_parameters', param.mklEps, 0, param.mklNorm );

% --- set kernels
nofKernels = size( Ks, 3 );
fprintf( 'setting %d kernels\n', nofKernels );
sg( 'set_kernel', 'COMBINED', param.cacheSize );
for( k = 1:nofKernels )
  K = Ks( idxTrn, idxTrn, k );
  sg( 'add_kernel', 1, 'CUSTOM', K, 'FULL' );
end;
fprintf( '\n' );



% === Prep Info Struct

info0 = [];
info0.idxTrn = idxTrn;
info0.idxVal = idxVal;
info0.idxTst = idxTst;

