function trSTp = compuTrSTp(K)
% Author: Xinwang Liu
n=size(K,1);
d=size(K,3);
trSTp = zeros(d,1);

for i=1:d
    trSTp(i) = (trace(K(:,:,i)) - sum(sum(K(:,:,i)))/n);
end
trSTp = trSTp(:)';