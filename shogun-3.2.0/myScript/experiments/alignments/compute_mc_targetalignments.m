function talignments = compute_mc_targetalignments(K,y);

% init
talignments = zeros(1,size(K,3));

y = y(:);
y = y-1;
ind = find(y==2);
y(ind)= complex(0,1);
ind = find(y==0);
y(ind)= complex(0,-1);

Y = y*y';

kcenter(Y);
kcenter(K);


% for each pair of kernels
for i=1:size(K,3)
  talignments(i) = sum(sum(K(:,:,i).*Y)) / sqrt( sum(sum(K(:,:,i).*K(:,:,i))) * sum(sum(Y.*Y)) );
  %talignments(i) = dotprod_spec(K(:,:,i),Y) / sqrt(dotprod_spec(K(:,:,i),K(:,:,i))*dotprod_spec(Y,Y));
end

