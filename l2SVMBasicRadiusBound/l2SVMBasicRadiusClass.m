function [radius,radiusp,xsup,w,d,pos,alpha,obj]=l2SVMBasicRadiusClass(x,y,C,lambda0, K0, Sigma, verbose,alphainit)

n=size(K0,1);	
[radius, radiusp]= basicRadius(K0,Sigma,lambda0,verbose);

K=sumKbeta(K0,Sigma);
%----------------------------------------------------------------------
%      monqp(H,b,c) solves the quadratic programming problem:
% 
%    min 0.5*x'Hx - d'x   subject to:  A'x = b  and  0 <= x <= c 
%     x    
%----------------------------------------------------------------------
Km = (K + eye(size(K))*radius/C);
H=(y*y').*Km;

f = ones(n,1);
A=y;
b=0;

if verbose ~= 0 
    disp('in QP'); 
end;

[alpha , d, pos] = monqpCinfty(H,f,A,b,lambda0,verbose,x,Km,alphainit);  

if verbose ~= 0 
    disp('out QP'); 
end;

alphaall=zeros(n,1);
alphaall(pos)=alpha;
obj=-0.5*alphaall'*H*alphaall +f'*alphaall;

if ~isempty(x)
    xsup=x(pos,:);
else
    xsup=[];
end
ysup=y(pos);
w=(alpha.*ysup);

