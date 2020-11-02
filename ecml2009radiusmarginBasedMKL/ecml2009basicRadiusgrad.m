function [grad] = ecml2009basicRadiusgrad(K,indsup,Alpsup,C,radiusp)

d=size(K,3);
grad=zeros(d,1);
for k=1:d  
    grad(k) = -0.5* Alpsup'* K(indsup,indsup,k)* Alpsup - 0.5*(radiusp(k)/C)*(Alpsup'*Alpsup); 
end;
