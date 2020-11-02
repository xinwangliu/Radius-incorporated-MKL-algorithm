function [grad] = l2SVMBasicRadiusMarginBoundGrad(K,indsup,Alpsup,C,radiusp)
%Author : Xinwang Liu
d=size(K,3);
n=size(K,1);
Kt=eye(n);

grad=zeros(d,1)';

for k=1:d;
    %  grad(k) = - 0.5*Alpsup'*Kaux(indsup,indsup)*(Alpsup)  ;
    grad(k) = - 0.5 * Alpsup'*( K(indsup,indsup,k) + Kt(indsup,indsup)*radiusp(k)/C )* Alpsup; 
end