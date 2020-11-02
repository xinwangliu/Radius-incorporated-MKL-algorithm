function [grad] = l2SVMradiusMarginBoundGrad(K,indsup,Alpsup,obj,radius,optimalBeta)
%Author : Xinwang Liu
d=size(K,3);
grad=zeros(d,1)';

for k=1:d;
    %  grad(k) = - 0.5*Alpsup'*Kaux(indsup,indsup)*(Alpsup)  ;
    grad(k) = - 0.5 * radius * Alpsup'* K(indsup,indsup,k) * Alpsup...
        +  ( sum(optimalBeta.*diag(K(:,:,k))) -   sum(sum((optimalBeta*optimalBeta').*K(:,:,k))) )*obj;   
end