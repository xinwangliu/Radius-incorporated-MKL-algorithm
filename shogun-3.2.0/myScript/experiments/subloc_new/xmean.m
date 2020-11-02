
function [ y ] = xmean( x, d );

if( ~exist('d','var') )
  d = 1;
end;

if( ~any( isnan(x(:)) ) )
  y = mean( x, d );
  return;
end;

x0 = x;
x0( isnan(x) ) = 0;
sums = sum( x0, d );
counts = sum( ~isnan(x), d );
ASSERT( sums(counts==0) == 0 );
counts( counts == 0 ) = 1;
y = sums ./ counts;

