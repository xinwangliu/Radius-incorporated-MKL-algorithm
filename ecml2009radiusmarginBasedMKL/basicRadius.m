function [radiusp]= basicRadius(K)

d=size(K,3);
n=size(K,1);
lambda0 = 1e-8;
radiusp = zeros(d,1);
for p = 1:d
    K(:,:,p) = (K(:,:,p)+K(:,:,p)')/2 + lambda0*eye(n);
    H= 2 * K(:,:,p);
    f= diag(K(:,:,p));
    A=ones(n,1);
    b=1;
    B=1;
    %% [xnew, lambda, pos,mu] = mymonqp(H,c,A,b,C,l)

    [beta, lambdat, pos] = mymonqp(H,f,A,b,B);
    
    optimalBeta=zeros(n,1);
    optimalBeta(pos)=beta;
    
    radiusp(p) = -0.5*optimalBeta'*H*optimalBeta + f'*optimalBeta;
end