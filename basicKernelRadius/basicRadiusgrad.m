function [grad] = basicRadiusgrad(K,indsup,Alpsup,radius,radiusp,Sigma)


d=size(K,3);
sumK=sumKbeta(K,Sigma);

grad=zeros(d,1)';
for k=1:d  
    grad(k) = - Alpsup'*K(indsup,indsup,k)*(Alpsup)/(2*radius)...
        + (Alpsup'*sumK(indsup,indsup)*Alpsup)*radiusp(k)/(2*radius^2);   
end;