function [grad] = l2trStGrad(K,indsup,Alpsup)


d=size(K,3);

grad=zeros(d,1)';
for k=1:d
     grad(k) = - 0.5*Alpsup'*K(indsup,indsup,k)*(Alpsup)  ;  
end;