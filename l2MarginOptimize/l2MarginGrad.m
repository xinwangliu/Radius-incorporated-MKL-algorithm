function [grad] = l2MarginGrad(K,indsup,Alpsup,Sigma)
%Author : Xinwang Liu

d=size(K,3);
grad(1)= + sum(Alpsup.*Alpsup) * Sigma(1)^2;
for k=2:d+1    
    grad(k) = - Alpsup'* K(indsup,indsup,k-1) * Alpsup ;
end

