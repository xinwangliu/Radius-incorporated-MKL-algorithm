function [KF,Sigma,Alpsup,w0,pos,obj] = myabsentmultikernelbinaryclassification(K,y,C,S,qnorm)

num = size(K,1);
nbkernel = size(K,3);
alpha0 = 1e-6;
%% S: num0*m, each column indicates the indices of absent samples
%% initialize kernel weights
Sigma = (1/nbkernel)^(1/qnorm)*ones(nbkernel,1);
%% initialize base kernels with zeros
KF = feval('algorithm2',K,S);
%% combining the base kernels
Kmatrix  = mycombFun(KF,Sigma);

loop = 1;
iter = 1;
while loop
    %%Update alpha with given Sigma & K
    %% [Alpsup,w0,pos,obj0] = mySVMclass(y,C,Kmatrix);
    [Alpsup,w0,pos,obj(iter)] = mySVMclass(y,C,Kmatrix);
    %-----------------------------------------
    % Update Sigma
    [grad] = gradsvmclass_sparse(KF,pos,Alpsup);
    norm_ker2 = -2*grad;
    norm_ker = sqrt(norm_ker2);
    norm_ker = Sigma.* norm_ker;
    % calculate Sigma
    Sigmaold = Sigma;
    % note!!!
    %
    if qnorm == 1
        Sigma = norm_ker/sum(norm_ker);
    else
        q = qnorm;
        norm_ker_q = norm_ker.^(2/(1+q));
        Sigma = norm_ker_q/sum(norm_ker_q.^q).^(1/q);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
   %% updata base kernels
    KF = zeros(num,num,nbkernel);
    for p =1:nbkernel
        if isempty(S{p}.indx)
            KF(:,:,p) = K(:,:,p);
        else
            %%[Kx] = myCalculateKx(Alpsup,pos,nbsv,num);
            Kx = myCalculateBinaryKx(Alpsup,pos,num);
            mis_indx = S{p}.indx;
            obs_indx = setdiff(1:num,mis_indx);
            KF(:,:,p) = absentKernelImputation(Kx,K(obs_indx,obs_indx,p),mis_indx,alpha0);
        end
    end
    Kmatrix = sumKbeta(KF,Sigma);
    %----------------------------------------------------
    % check variation of Sigma conditions
    %----------------------------------------------------
    if  max(abs(Sigma - Sigmaold))<1e-3 || ( iter>2 && (abs((obj(iter-1)-obj(iter))/obj(iter))<1e-4...
            ||iter>=100) )
        loop = 0;
        fprintf(1,'variation convergence criteria reached \n');
    end
    iter = iter + 1;
end