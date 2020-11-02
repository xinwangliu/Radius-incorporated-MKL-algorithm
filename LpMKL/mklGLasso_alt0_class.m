function [Sigma,Alpsup,w0,pos,obj] = mklGLasso_alt0_class(K,y,C,qnorm,verbose)

% USAGE [Sigma,Alpsup,w0,pos,history,obj,status] = mklsvm(K,yapp,C,option,verbose)
%
% Inputs
%
% K         : NxNxD matrix containing all the Gram Matrix
% C         : SVM hyperparameter
% option    : mkl algorithm hyperparameter
%
%       option.nbitermax : maximal number of iterations (default 1000)
%       option.algo      : selecting algorithm svmclass (default) or svmreg
%       option.seuil     : threshold for zeroing kernel coefficients
%                          (default 1e-12) during the first
%                          option.seuilitermax
%       option.seuilitermax : threshold weigth when the nb of iteration is
%                           lower than this value
%       option.Sigmainit  : initial kernel coefficient (default average)
%       option.alphainit : initial Lagrangian coefficient
%       option.lambdareg : ridge regularization added to kernel's diagonal
%                          for SVM (default 1e-10)
%       option.verbosesvm : verbosity of SVM (default 0) see svmclass or
%                           svmreg
%       option.svmreg_epsilon : epsilon for SVM regression (mandatory for
%                         svm regression)
%       option.numericalprecision    : force to 0 weigth lower than this
%                                     value (default = eps)
%
% Outputs
%
% Sigma         : the weigths
% Alpsup        : the weigthed lagrangian of the support vectors
% w0            : the bias
% pos           : the indices of SV
% history        : history of the weigths
% obj           : objective value
% status        : output status (sucessful or max iter)


nbkernel=size(K,3);
num = size(K,1);
Sigma=ones(nbkernel,1)/nbkernel;
%------------------------------------------------------------------------------%
% Initialize
%------------------------------------------------------------------------------%
nloop = 0;
loop = 1;
nloopmax = 100;
%-----------------------------------------
% Initializing SVM
%------------------------------------------
Kmatrix = sumKbeta(K,Sigma.^qnorm);
% initial a solution for future use
[Alpsup,w0,pos,obj0] = mySVMclass(y,C,Kmatrix);
while loop
    
    nloop = nloop+1;
    %-----------------------------------------
    % Update weigths Sigma
    %----------------------------------------
    %% [grad] = gradsvmclass_sparse(gamma,K,indsup,Alpsup)
    [grad] = gradsvmclass_sparse(K,pos,Alpsup);
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
    
    %     %----------------------------------------------------
    %     % Run SVM
    %     %----------------------------------------------------
    Kmatrix=sumKbeta(K,Sigma.^qnorm);
    [Alpsup,w0,pos,obj(nloop)] = mySVMclass(y,C,Kmatrix);

%     % calculate duality gap
%     normek=-grad ; % 0.5*Alpsup'*K(pos,pos,i)*Alpsup;
%     primalobj= C*sum(max( 1-y.*(Kmatrix(:,pos)*Alpsup + w0),0));
%     %posnz = find(Sigma > option.threshold);
%     %primalobj = primalobj+0.5*sum(norm_ker2(posnz).*Sigma(posnz));
%     primalobj = primalobj + sum(normek); %0.5*norm_ker2*Sigma';
%     if qnorm==1
%         dualitygap=(obj(nloop) + max(normek) - sum(abs(Alpsup)))/obj(nloop);       
%     else   
%         dualitygap= (  primalobj - obj(nloop) ) /primalobj;      
%     end
    %------------------------------------------
    %  verbosity
    %------------------------------------------
    %
    if verbose
        if nloop == 1 || rem(nloop,10)==0
            fprintf('--------------------------------------------------\n');
            fprintf('Iter | Obj.    | DiffSigmaSs  |\n');
            fprintf('--------------------------------------------------\n');
        end;
        fprintf('%d   | %8.4f | %6.4f  \n',[nloop obj(nloop)   max(abs(Sigma-Sigmaold)) ]);
    end
    %----------------------------------------------------
    % check variation of Sigma conditions
    %----------------------------------------------------
    if  max(abs(Sigma - Sigmaold))<1e-4 || ( nloop>2 && abs((obj(nloop-1)-obj(nloop))/obj(nloop))<1e-4 )
        loop = 0;
        fprintf(1,'variation convergence criteria reached \n');
    end
    %-----------------------------------------------------
    % check nbiteration conditions
    %----------------------------------------------------
    if nloop>=nloopmax
        loop = 0;
        fprintf(1,'maximum number of iterations reached\n')
    end
end