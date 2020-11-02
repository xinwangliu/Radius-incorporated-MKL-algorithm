function talignments = compute_targetalignments(K,y,Y);

% init
talignments = zeros(1,size(K,3));

if max(max(size(y)))>1
  if size(y,2)>size(y,1), y=y'; end
  Y = y*y';
end

kcenter(Y);
kcenter(K);

% for each pair of kernels
for i=1:size(K,3)
  talignments(i) = sum(sum(K(:,:,i).*Y)) / sqrt( sum(sum(K(:,:,i).*K(:,:,i))) * sum(sum(Y.*Y)) );
  %talignments(i) = dotprod_spec(K(:,:,i),Y) / sqrt(dotprod_spec(K(:,:,i),K(:,:,i))*dotprod_spec(Y,Y));
end

