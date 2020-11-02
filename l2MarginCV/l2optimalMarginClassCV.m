function [xsup,w,d,pos,alpha,margin]=l2optimalMarginClassCV(x,y,C,lambda, kernel, kerneloption,verbose,alphainit)

K=svmkernel(x,kernel,kerneloption);
n=size(K,1);	

%----------------------------------------------------------------------
%      monqp(H,b,c) solves the quadratic programming problem:
% 
%    min 0.5*x'Hx - d'x   subject to:  A'x = b  and  0 <= x <= c 
%     x    
%----------------------------------------------------------------------
Km = (K + eye(size(K))/C);
H=(y*y').*Km;

f = ones(n,1);
A=y;
b=0;

if verbose ~= 0 
    disp('in QP'); 
end;

[alpha , d, pos] = monqpCinfty(H,f,A,b,lambda,verbose,x,Km,alphainit);  

if verbose ~= 0 
    disp('out QP'); 
end;

alphaall=zeros(n,1);
alphaall(pos)=alpha;

margin = alphaall'*H*alphaall;

if ~isempty(x)
    xsup=x(pos,:);
else
    xsup=[];
end
ysup=y(pos);
w=(alpha.*ysup);

