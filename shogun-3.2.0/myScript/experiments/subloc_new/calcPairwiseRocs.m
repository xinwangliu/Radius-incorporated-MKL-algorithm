
function [ aucs ] = calcPairwiseRocs( y, outs );

[ N, C ] = size( outs );
ASSERT( size(y) == [N,1] );
classes = unique( y );
ASSERT( length(classes) == C );

aucs = repmat( nan, 1, C*(C-1)/2 );
l = 0;
for( i1 = 1:C )
  y1 = classes( i1 );
  idx1 = find( y == y1 );
  for( i2 = (i1+1):C )
    y2 = classes( i2 );
    idx2 = find( y == y2 );
    idx = [ idx1 ; idx2 ];
    out12 = outs(idx,i1) - outs(idx,i2);
    y12 = (y(idx)==y1) - (y(idx)==y2);
    res = calcRoc( y12, out12 );
    l = l + 1;
    aucs(l) = res.aucRoc;
  end;
end;

