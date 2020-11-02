function [grad] = l2RadiusMarginBoundGradCV(K,indsup,Alpsup,margin,radius,optimalBeta)
%Author : Xinwang Liu
d=size(K,3);
grad=zeros(d,1)';

for k=1:d    
    grad(k) =  - Alpsup'* K(indsup,indsup,k) * Alpsup * radius...
        +  ( sum(optimalBeta.*diag(K(:,:,k))) - sum(sum((optimalBeta*optimalBeta').*K(:,:,k))) ) * margin;   
end