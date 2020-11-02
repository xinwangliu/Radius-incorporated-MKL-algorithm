
function [ res ] = collectResults( datasetName, dirName, measure, splits );

%%% datasetName = 'plant';
%%% datasetName = 'nonpl';
%%% measure = 'MCC';
%%% measure = 'AUC';

resDir = [ 'res/' dirName ];
dataDir = [ 'data/' datasetName '/' ];

labelFileName = [ dataDir 'label_' datasetName '.mat' ];
load( labelFileName, 'y' );

if( ~exist('splits','var') )
  splits = [];
end;
if( isempty(splits) )
  permsFileName = [ dataDir 'perms,' datasetName '.mat' ];
  load( permsFileName, 'perms' );
  [ nofPerms, N ] = size( perms );
  splits = 1:nofPerms;
end;

res = [];
l = 0;
%for( numPerm = 1:nofPerms )
for( numPerm = splits )
  resFileName = sprintf( '%s%s,split_%02d.mat', resDir, datasetName, numPerm );
  %fprintf( '%s\n', resFileName );
  if( exist(resFileName,'file') )
    load( resFileName, 'RES', 'INFO', 'param0' );
    if( any(isnan(RES(:))) )
      %continue;
    end;
    l = l + 1;
    switch( measure )
     case 'MCC'
      r = mean( RES, 3 );
     case 'AUC'
      idxTst = INFO{1,1}.idxTst;
      yTst = y( idxTst );
      [ nofNorms, nofCs, nofClasses ] = size( RES );
      r = repmat( nan, nofNorms, nofCs );
      for( i = 1:nofNorms )
	for( j = 1:nofCs )
	  info = INFO{ i, j };
	  if( isempty(info) )
	    continue;
	  end;
	  ASSERT( info.idxTst == idxTst );
	  outTst = info.outs( idxTst, : );
	  aucs = calcPairwiseRocs( yTst, outTst );
	  r( i, j ) = mean( aucs );
	end;
      end;
     otherwise
      error( '???' );
    end;
    res(:,:,l) = r;
    fprintf( 'split %02d:  %3d values ready\n', numPerm, sum(isfinite(r(:))) );
  end;
end;



% === results

% --- norm strs
[ nofNorms, nofCs, nofClasses ] = size( RES );
xLabStrs = {};
for( i = 1:nofNorms )
  xLabStrs{i} = num2str( param0.mklNorms(i) );
end;

% --- table header
fprintf( '\n' );
fprintf( '{\\bf norm}' );
fprintf( ' & %s', xLabStrs{:} );
fprintf( ' \\\\\n' );

% --- different evaluation modes
for( i = 2 )

  % --- eval
  S = [];
  switch( i )
   case 1
    t = 'maxAvg';
    T = max( xmean(res,3), [], 2 );
   case 2
    t = 'avgMax';
    T = xmean( max(res,[],2), 3 );
    S = std( squeeze( max(res,[],2) )' );
   otherwise
    error( '???' );
  end;
  ASSERT( 0 <= T | isnan(T) );
  ASSERT( T <= 1 | isnan(T) );

  % --- figure
  figure;
  bar( 1-T );
  title( datasetName );
  xlabel( 'MKL norm' );
  ylabel( [ '1 - ' t measure ] );
  set( gca, 'XTickLabel', xLabStrs );
  plotFileName = [ resDir datasetName ',' t measure '.ps' ];
  print( '-dpsc', plotFileName );

  % --- table
  fprintf( '\\hline {\\bf 1 - %s %s [\\%%]}', t, measure );
  fprintf( ' & $%.2f$', 100*(1-T) );
  fprintf( ' \\\\\n' );
  if( ~isempty(S) )
    fprintf( '\\hline {\\bf SD(%s) [\\%%]}', measure );
    fprintf( ' & $%.2f$', 100*S );
    fprintf( ' \\\\\n' );
  end;
  
end;

fprintf( '\n' );


