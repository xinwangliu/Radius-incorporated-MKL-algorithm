
function [ res ] = calcRoc( y, out, opt );

% === init
n = length( out );
nofPos = sum( y == +1 );
nofNeg = sum( y == -1 );
% --- check input
ASSERT( size(out) == size(y) );
ASSERT( nofPos+nofNeg == n );
ASSERT( nofPos > 0 );
ASSERT( nofNeg > 0 );
ASSERT( ~isnan(out) );
% --- default options
if( ~ exist('opt','var') )
  opt = [];
end;
opt = setDefault( opt, 'b', 0.0 );
opt = setDefault( opt, 'plotRoc', 0 );
opt = setDefault( opt, 'plotFdr', 0 );
opt = setDefault( opt, 'title', '' );
opt = setDefault( opt, 'grid', 1 );
opt = setDefault( opt, 'granularity', 1000 );
% --- split into pos and neg
outPos = out( y == +1 );
outNeg = out( y == -1 );
granularity = opt.granularity;

% === compute ROC
if( 0 )
[ dummy, perm ] = sort( out );
clear dummy;
TP = nofPos;
TN = 0;
FP = nofNeg;
FN = 0;
TPRs = nan*zeros( n+1, 1 );
FPRs = nan*zeros( n+1, 1 );
PPVs = nan*zeros( n+1, 1 );
TPRs(1) = 1;
FPRs(1) = 1;
PPVs(1) = 0;
for( i = 1:n )
  j = perm( i );
  if( y(j) == -1 )
    FP = FP - 1;
    TN = TN + 1;
  else
    TP = TP - 1;
    FN = FN + 1;
  end;
  sn = TP / ( TP + FN );
  sp = TN / ( TN + FP );
  TPRs(i+1) = sn;
  FPRs(i+1) = 1 - sp;
  if( TP+FP == 0 )
    PPVs(i+1) = 1;
  else
    PPVs(i+1) = TP / ( TP + FP );
  end;
end;
ASSERT( TP == 0 );
ASSERT( TN == nofNeg );
ASSERT( FP == 0 );
ASSERT( FN == nofPos );
end;

% === compute ROC
thresholds = -sort( -outPos );
if( nofPos < 2*granularity )
  p = nofPos;
else
  p = granularity;
end;
TPRs = (1:p) / p;
FPRs = nan*zeros( p, 1 );
PPVs = nan*zeros( p, 1 );
for( i = 1:p )
  tpr = TPRs( i );
  threshold = thresholds( ceil(nofPos/p*i) );
  tp = sum( outPos >= threshold );
  fp = sum( outNeg >= threshold );
  ASSERT( tp/nofPos >= tpr );
  % --- correct for ties
  nofPosMatch = sum( outPos == threshold );
  nofNegMatch = sum( outNeg == threshold );
  TP = nofPos * tpr;
  gapPos = round( tp - TP );
  ASSERT( gapPos >= 0 );
  if( gapPos > 0 )
    frac = gapPos / nofPosMatch;
    gapNeg = frac * nofNegMatch;
    tp = TP;
    fp = fp - gapNeg;
  end;
  % --- store point
  FPRs(i) = fp / nofNeg;
  PPVs(i) = tp / ( tp + fp );
end;

% === evaluate for given bias
TP = sum( y==+1 & out>opt.b );
TN = sum( y==-1 & out<opt.b );
FN = sum( y==+1 & out<opt.b );
FP = sum( y==-1 & out>opt.b );
res.sn = TP / ( TP + FN );
res.sp = TN / ( TN + FP );
if( TP+FP == 0 )
  res.ppv = 0;
else
  res.ppv = TP / ( TP + FP );
end;

% === return results
res.TPRs = TPRs;
res.FPRs = FPRs;
res.PPVs = PPVs;
res.aucRoc = 1 - mean( FPRs );
res.aucPpv = mean( PPVs );

% === plot ROC
if( opt.plotRoc )
  res.rocPlot = figure;
  hold on;
  plot( FPRs, TPRs, 'LineWidth', 2 );
  plot( 1-res.sp, res.sn, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 10, 'Color', 'r', 'LineStyle', 'none' );
  line( [ 0 1 ], [ 0 1 ], 'Color', 'k' );
  if( ~ isempty( opt.title ) )
    title( sprintf( 'ROC, %s', opt.title ) );
  end;
  xlabel( 'false positive rate' );
  ylabel( 'true positive rate' );
  legendStr1 = sprintf( 'ROC, AUC=%.2f%%', 100*res.aucRoc );
  legendStr2 = sprintf( 'Sn=%.2f%%, Sp=%.2f%%', 100*res.sn, 100*res.sp );
  legendStrs = { legendStr1, legendStr2 };
  %legend( legendStrs, 'Location', 'SouthEast' );
  if( opt.grid )
    grid on;
  end;
end;

% === plot PPV vs SN
if( opt.plotFdr )
  res.fdrPlot = figure;
  hold on;
  plot( PPVs, TPRs, 'LineWidth', 2 );
  plot( res.ppv, res.sn, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 10, 'Color', 'r', 'LineStyle', 'none' );
  line( [ 0 1 ], [ 0 1 ], 'Color', 'k' );
  if( ~ isempty( opt.title ) )
    title( sprintf( 'PPV, %s', opt.title ) );
  end;
  xlabel( 'positive predictive value [1-FDR]' );
  ylabel( 'sensitivity (TPR)' );
  legendStr1 = sprintf( 'PPV, AUC=%.2f%%', 100*res.aucPpv );
  legendStr2 = sprintf( 'TPR=%.2f%%, PPV=%.2f%%', 100*res.sn, 100*res.ppv );
  legendStrs = { legendStr1, legendStr2 };
  %legend( legendStrs, 'Location', 'SouthEast' );
  if( opt.grid )
    grid on;
  end;
end;

