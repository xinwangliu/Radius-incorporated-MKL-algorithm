function [grad] = radiusMarginBoundGrad(K,Sigma,indsup,Alpsup,radius,optimalBeta)
%Author : Xinwang Liu
d=size(K,3);
sumK = sumKbeta(K,Sigma);
grad = zeros(1,d);
for k=1:d    
    grad(k) = - Alpsup'*K(indsup,indsup,k)*(Alpsup)/(2*radius)...
        + (Alpsup'*sumK(indsup,indsup)*Alpsup)*( sum(optimalBeta.*diag(K(:,:,k)))...
        - sum(sum((optimalBeta*optimalBeta').*K(:,:,k))))/(2*radius^2);  
end
