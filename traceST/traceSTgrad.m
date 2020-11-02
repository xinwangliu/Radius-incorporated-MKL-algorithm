function [grad] = traceSTgrad(K,indsup,Alpsup,trST,trSTp,Sigma)


d=size(K,3);
sumK=sumKbeta(K,Sigma);

grad=zeros(d,1)';
for k=1:d;
    %  grad(k) = - 0.5*Alpsup'*Kaux(indsup,indsup)*(Alpsup)  ;
    grad(k) = - Alpsup'*K(indsup,indsup,k)*(Alpsup)/(2*trST)...
        + (Alpsup'*sumK(indsup,indsup)*Alpsup)*trSTp(k)/(2*trST^2);   
end;