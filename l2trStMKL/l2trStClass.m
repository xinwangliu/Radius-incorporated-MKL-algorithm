function [xsup,w,d,pos,alpha,obj]=l2trStClass(x,y,lambda, K, Sigma,verbose,alphainit)

K0 = sumKbeta(K,Sigma);
n=size(K0,1);	

%----------------------------------------------------------------------
%      monqp(H,b) solves the quadratic programming problem:
% 
%    min 0.5*x'Hx - d'x   subject to:  A'x = b  and  0 <= x  
%     x    
%----------------------------------------------------------------------
H=(y*y').*K0;
f = ones(n,1);
A=y;
b=0;

if verbose ~= 0 
    disp('in QP'); 
end;

[alpha, d, pos] = monqpCinfty(H,f,A,b,lambda,verbose,x,K0,alphainit);  

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

