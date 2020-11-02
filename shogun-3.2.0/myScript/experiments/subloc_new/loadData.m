
function [ y, Ks, perms, info ] = loadData( datasetName );

dirName = [ 'data/' datasetName '/' ];
labelFileName = [ dirName 'label_' datasetName '.mat' ];
klistFileName = [ dirName 'klist_' datasetName '.mat' ];
permsFileName = [ dirName 'perms,' datasetName '.mat' ];

load( labelFileName, 'y' );
load( klistFileName, 'klist' );
load( permsFileName, 'perms' );

N = length( y );
nofKernels = length( klist );
ASSERT( size(perms,2) == N );

Ks = repmat( nan, [ N, N, nofKernels ] );
for( k = 1:nofKernels )
  kernFileName = [ dirName klist{k} ];
  load( kernFileName, 'K' );
  s = 1 / ( mean(diag(K)) - mean(K(:)) );
  K = s * K;
  Ks( :, :, k ) = K;
  scales(k) = s;
end;

info = [];
info.klist = klist;
info.scales = scales;

