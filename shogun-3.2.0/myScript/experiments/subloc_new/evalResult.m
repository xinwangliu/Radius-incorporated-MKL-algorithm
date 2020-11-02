

% === init

method = [ methodName '_' epsSetting ];
fprintf( '=== %s, %s ===\n', datasetName, method );

betaFileName = sprintf( 'res/BETA/%s,%s.mat', datasetName, method );
evalFileName = sprintf( 'res/EVAL/%s,%s.mat', datasetName, method );



% === collect results

if( ~exist(evalFileName,'file') )

Bs = [];  % betas
Gs = [];  % max gap    
Rs = [];  % max rel gap
Ts = [];  % train times

for( split = 1:nofSplits )

  resFileName = sprintf( 'res/%s/%s/%s,split_%02d.mat', datasetName, method, datasetName, split );
  fprintf( 'loading %s\n', resFileName );
  load( resFileName, 'INFO', 'RES', 'param0' );
  nofNorms = length( param0.mklNorms );
  nofCs = length( param0.logCs );
  [ nofKernels, nofClasses ] = size( INFO{1,1}.betas );

  B = repmat( nan, [ nofNorms, nofCs, nofKernels, nofClasses ] );
  G = repmat( nan, nofNorms, nofCs );
  R = repmat( nan, nofNorms, nofCs );
  T = repmat( nan, nofNorms, nofCs );
  
  for( iNorm = 1:nofNorms )
    mklNorm = param0.mklNorms( iNorm );
    for( iC = 1:nofCs )
      C = 2 ^ param0.logCs( iC );
      info = INFO{ iNorm, iC };
      gaps = info.objs(1,:) - info.objs(2,:);
      if( any(gaps<=0) )
	fprintf( 'WARNING: duality gap %e\n', min(gaps) );
      end;
      %ASSERT( gaps >= 0 );
      relGaps = gaps ./ info.objs(1,:);
      [ dummy, k ] = max( abs(relGaps) );
      maxRelGap = relGaps(k);
      B( iNorm, iC, :, : ) = info.betas;
      G( iNorm, iC ) = mean( gaps );
      R( iNorm, iC ) = maxRelGap;
      T( iNorm, iC ) = sum( info.trainTimes );
    end;
  end;

  Bs(:,:,:,:,split) = B;
  Gs(:,:,split) = G;
  Rs(:,:,split) = R;
  Ts(:,:,split) = T;
end;
fprintf( '\n' );

save( betaFileName, 'param0', 'Bs' );
save( evalFileName, 'param0', 'Gs', 'Rs', 'Ts' );

end;



% === plot results

load( evalFileName, 'param0', 'Gs', 'Rs', 'Ts' );
nofNorms = length( param0.mklNorms );
nofCs = length( param0.logCs );


% --- mean absolute duality gaps

figure;
imagesc( mean(Gs,3) );
colorbar( 'vert' );
titleStr = sprintf( 'avg. abs. duality gaps -- %s -- %s', datasetName, method );
title( titleStr );
set( gca, 'XTickLabel', param0.logCs );
set( gca, 'YTickLabel', param0.mklNorms );
xlabel( 'log_2(C)' );
ylabel( 'MKL norm' );
plotFileName = sprintf( 'res/PLOTS/avgAbsGaps,%s,%s.ps', datasetName, method );
print( '-dpsc', plotFileName );
close;


% --- max relative duality gaps

figure;
imagesc( max(Rs,[],3) );
colorbar( 'vert' );
titleStr = sprintf( 'max. rel. duality gaps -- %s -- %s', datasetName, method );
title( titleStr );
set( gca, 'XTickLabel', param0.logCs );
set( gca, 'YTickLabel', param0.mklNorms );
xlabel( 'log_2(C)' );
ylabel( 'MKL norm' );
plotFileName = sprintf( 'res/PLOTS/maxRelGaps,%s,%s.ps', datasetName, method );
print( '-dpsc', plotFileName );


% --- training times

avgTrainTime = mean(Ts(:));
T = mean( Ts, 3 );
tMax = max( T(:) );
tMin = min( T(:) );

axStrsMkl = {};
for( i = 1:nofNorms )
  axStrsMkl{i} = num2str( param0.mklNorms(i) );
end;

legStrsMkl = {};
for( i = 1:nofNorms )
  legStrsMkl{i} = [ 'p=', num2str( param0.mklNorms(i) ) ];
end;

legStrsC = {};
for( i = 1:nofCs )
  legStrsC{i} = sprintf( 'C=2^%d', param0.logCs(i) );
end;

% --  by C

if( 1 )
  figure;
  hold on;
  plot( param0.logCs, T', 'x-' );
  plot( param0.logCs, mean(T,1), 'Color', 'k', 'LineWidth', 3 );
  set( gca, 'YScale', 'log' );
  ax = axis;
  axis( [ ax(1) ax(2) tMin tMax ] );
  grid on;
  titleStr = sprintf( 'training times -- %s -- %s (avg %.0f secs)', datasetName, method, avgTrainTime );
  title( titleStr );
  xlabel( 'log_2(C)' );
  ylabel( 'training time' );
  legend( { legStrsMkl{:}, 'avg' }, 3 );
  plotFileName = sprintf( 'res/PLOTS/trainTimes,C,%s,%s.ps', datasetName, method );
  print( '-dpsc', plotFileName );
  close;
end;

% -- by p

figure;
hold on;
plot( T, 'x-' );
plot( mean(T,2), 'Color', 'k', 'LineWidth', 3 );
set( gca, 'YScale', 'log' );
ax = axis;
axis( [ ax(1) ax(2) tMin tMax ] );
grid on;
titleStr = sprintf( 'training times -- %s -- %s (avg %.0f secs)', datasetName, method, avgTrainTime );
title( titleStr );
set( gca, 'XTickLabel', axStrsMkl );
xlabel( 'MKL norm' );
ylabel( 'training time' );
legend( { legStrsC{:}, 'avg' }, 3 );
plotFileName = sprintf( 'res/PLOTS/trainTimes,p,%s,%s.ps', datasetName, method );
print( '-dpsc', plotFileName );

