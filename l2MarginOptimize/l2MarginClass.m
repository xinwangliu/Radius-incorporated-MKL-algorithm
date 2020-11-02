function [xsup,w,d,pos,alpha,obj] = l2MarginClass(x,y,lambda,K0,Sigma,verbose,alphainit)

weight = Sigma(2:end);
K = sumKbeta(K0,weight);

Km = K + eye(size(K))*Sigma(1);
n=size(Km,1);	

%----------------------------------------------------------------------
%      monqp(H,b,c) solves the quadratic programming problem:
% 
%    min 0.5*x'Hx - d'x   subject to:  A'x = b  and  0 <= x <= c 
%     x    
%----------------------------------------------------------------------
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
%margin = -0.5*alphaall'*H*alphaall +f'*alphaall;

obj = alphaall'*H*alphaall;

if ~isempty(x)
    xsup=x(pos,:);
else
    xsup=[];
end
ysup=y(pos);
w=(alpha.*ysup);
