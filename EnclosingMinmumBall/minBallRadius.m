function [radius, optimalBeta]= minBallRadius(K,lambda,verbose)

n=size(K,1);
H=2*K;
f=diag(K);
A=ones(n,1);
b=1;
C=ones(n,1);
L=zeros(n,1);

[beta, lambda, pos] = monqp(H,f,A,b,C,lambda,verbose);

optimalBeta=zeros(n,1);
optimalBeta(pos)=beta;

radius = -0.5*optimalBeta'*H*optimalBeta + f'*optimalBeta;

% 
