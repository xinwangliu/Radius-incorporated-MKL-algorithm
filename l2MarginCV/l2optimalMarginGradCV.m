function [grad] = l2optimalMarginGradCV(K,indsup,Alpsup)
%Author : Xinwang Liu
d=size(K,3);
grad=zeros(d,1)';

for k=1:d    
    grad(k) =  - Alpsup'* K(indsup,indsup,k) * Alpsup;  
end