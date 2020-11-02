function [trSTp]= traceST(K)

d=size(K,3);
n=size(K,1);

trSTp=zeros(d,1);
for p=1:d
    trSTp(p) = (trace(K(:,:,p)) - sum(sum(K(:,:,p)))/n)/n;
end
