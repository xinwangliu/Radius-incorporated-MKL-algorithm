function [grad] = gradMiniBallClass(K,indsup,Alpsup,radius,optimalBeta,Sigma)


d=size(K,3);
sumK=sumKbeta(K,Sigma);

grad=zeros(d,1)';
for k=1:d;
    %  grad(k) = - 0.5*Alpsup'*Kaux(indsup,indsup)*(Alpsup)  ;
    grad(k) = - Alpsup'*K(indsup,indsup,k)*(Alpsup)/(2*radius)...
        + (Alpsup'*sumK(indsup,indsup)*Alpsup)*( sum(optimalBeta.*diag(K(:,:,k)))...
        - sum(sum((optimalBeta*optimalBeta').*K(:,:,k))))/(2*radius^2);   
end;
