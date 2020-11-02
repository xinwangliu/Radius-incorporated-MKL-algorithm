
function [ RES, INFO ] = trainAndTestSeveral( param0 );


% === init
datasetName = param0.datasetName;
numPerm = param0.numPerm;
fprintf( '\n' );
fprintf( '=== %s, %s, %s, split %02d ===\n', datasetName, param0.optName, param0.epsSetting, numPerm );
fprintf( '\n' );
nofNorms = length( param0.mklNorms );
nofCs = length( param0.logCs );

% --- load data
[ y, Ks, perms, dataInfo ] = loadData( param0.datasetName );
classes = unique( y );
nofClasses = length( classes );

% --- res struct: already computed?
resFileName = sprintf( '%s%s,split_%02d.mat', param0.resDir, datasetName, numPerm );
param0.resFileName = resFileName;
if( ~exist(resFileName,'file') )
  % --- prep 
  RES = repmat( nan, [ nofNorms, nofCs, nofClasses ] );
  INFO = cell( nofNorms, nofCs );
  save( resFileName, 'RES', 'INFO', 'param0' );
else
  % --- load
  fprintf( 'loading %s\n', resFileName );
  load( resFileName, 'RES', 'INFO' );
end;

% === compute
[ info0 ] = prepForMkl( param0, y, Ks, perms );
param0.idxTrn = info0.idxTrn;
param0.idxTst = info0.idxTst;

% --- go thru models
for( iNorm = 1:nofNorms )
  mklNorm = param0.mklNorms( iNorm );
  for( iC = 1:nofCs )
    C = 2 ^ param0.logCs( iC );
    fprintf( '--- norm %5.2f, C %5.2f ---\n', mklNorm, C );
    load( resFileName, 'RES', 'INFO' );
    if( any( RES(iNorm,iC,:) < 0 ) )
      fprintf( 'previous computation ongoing or aborted.\n' );
      if( 1 )
	fprintf( 'starting recomputation.\n' );
	RES(iNorm,iC,:) = nan;
      else
	fprintf( 'skipping computation.\n' );
	fprintf( '\n' );
	continue;
      end;
    end;
    if( ~any(isnan( RES(iNorm,iC,:) )) )
      fprintf( 'already computed - skipping.\n' );
      fprintf( '\n' );
      continue;
    end;
    RES( iNorm, iC, : ) = -1;
    save( resFileName, 'RES', 'INFO', 'param0' );
    
    % --- settings
    param = param0;
    param.svmC = C;
    param.mklNorm = mklNorm;
    if( mklNorm == 1 && length(param0.optName) >= 6 && strcmp(param0.optName(1:6),'NEWTON') )
      param.mklSolver = 'CPLEX';
    end;
    sg( 'c', param.svmC );
    sg( 'set_solver', param.mklSolver );
    sg( 'mkl_parameters', param.mklEps, 0, param.mklNorm );
    %%% if( mklNorm == inf )
    %%%   fprintf( 'fixing betas to results of %.3f-norm\n', param0.mklNorms(iNorm-1) );
    %%%   param.fixBetas = INFO{ iNorm-1, iC }.betas;
    %%% end;
    % --- assess model
    [ res, info ] = trainAndTest( param, y, Ks, perms, info0 );
    load( resFileName, 'RES', 'INFO' );
    if( all( RES(iNorm,iC,:) == -1 ) )
      RES(  iNorm, iC, : ) = res.MCC;
      INFO{ iNorm, iC } = info;
      save( resFileName, 'RES', 'INFO', 'param0' );
    else
      fprintf( 'WARNING: computation clash -- reporting inserted results\n' );
      res = [];
      res.MCC = RES( iNorm, iC, : );
      info = INFO{ iNorm, iC };
    end;
    gaps = info.objs(1,:) - info.objs(2,:);
    relGaps = gaps ./ info.objs(1,:);
    fprintf( 'max gap:  %.2f / %.1f => %.3f rel \n', max(gaps), mean(info.objs(1,:)), max(relGaps) );
    fprintf( 'avg MCC:  %5.2f%% \n', 100*mean(res.MCC) );
    fprintf( '\n' );
  end;
end;

% --- finish
save( resFileName, 'RES', 'INFO', 'param0' );
fprintf( '\n' );
fprintf( '=== DONE: %s, split %02d ===\n', datasetName, numPerm );
fprintf( '\n' );

