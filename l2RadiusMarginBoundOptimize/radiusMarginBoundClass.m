function [radius, optimalBeta,xsup,w,d,pos,alpha,margin,obj] = radiusMarginBoundClass(x,y,lambda,K,Sigma,verbose,alphainit) 

% compute radius;
[radius, optimalBeta]= radiusCompu(K,Sigma,lambda,verbose);

Km = sumKbeta(K,Sigma);
n=size(Km,1);	

%----------------------------------------------------------------------
%      monqp(H,b,c) solves the quadratic programming problem:
% 
%    min 0.5*x'Hx - d'x   subject to:  A'x = b  and  0 <= x <= c 
%     x    
%----------------------------------------------------------------------
H=(y*y').*Km/radius;
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

Km0 = sumKbeta(K,Sigma(1:end-1));
margin = (alphaall'*Km0*alphaall)/(radius^2);

if ~isempty(x)
    xsup=x(pos,:);
else
    xsup=[];
end
ysup=y(pos);
w=(alpha.*ysup);

obj = radius *margin;
