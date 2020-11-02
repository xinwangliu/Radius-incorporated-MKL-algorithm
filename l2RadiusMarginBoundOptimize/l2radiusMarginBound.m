function [radius,Sigma,Alpsup,w0,pos,obj,status] = l2radiusMarginBound(K,yapp,option,verbose)
%% Author: Xinwang Liu
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
%       option.sigmainit  : initial kernel coefficient (default average)
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

if ~isempty(K)
    if size(K,3)>1
        nbkernel=size(K,3);
        if option.efficientkernel==1
            K = build_efficientK(K);
        end;
    elseif option.efficientkernel==1 & isstruct(K);
        nbkernel=K.nbkernel;     
    end;
else
    error('No kernels defined ...');
end;

if ~isfield(option,'nbitermax');
    nloopmax=1000;
else
    nloopmax=option.nbitermax;
end;
if ~isfield(option,'algo');
    option.algo='svmclass';
end;
if ~isfield(option,'seuil');
    seuil=1e-12;
else
    seuil=option.seuil;
end
if ~isfield(option,'seuilitermax')
    option.seuilitermax=20;
end;
if ~isfield(option,'seuildiffsigma');
    option.seuildiffsigma=1e-5;
end
if ~isfield(option,'seuildiffconstraint');
    option.seuildiffconstraint=0.05;
end

if ~isfield(option,'lambdareg');
    lambdareg=1e-10;
    option.lambdareg=1e-10;
else
    lambdareg=option.lambdareg;
end


if ~isfield(option,'numericalprecision');
    option.numericalprecision=0;
end;

if ~isfield(option,'verbosesvm');
    verbosesvm=0;
    option.verbosesvm=0;
else
    verbosesvm=option.verbosesvm;
end

if ~isfield(option,'sigmainit');
    Sigma = ones(1,nbkernel)/nbkernel;
else
    Sigma=option.sigmainit ;
    ind=find(Sigma==0);
end;


if isfield(option,'alphainit');
    alphainit=option.alphainit;
else
    alphainit=[];
end;




%--------------------------------------------------------------------------------
% Options used in subroutines
%--------------------------------------------------------------------------------
if ~isfield(option,'goldensearch_deltmax');
    option.goldensearch_deltmax=1e-1;
end
if ~isfield(option,'goldensearchmax');
    optiongoldensearchmax=1e-8;
end;
if ~isfield(option,'firstbasevariable');
    option.firstbasevariable='first';
end;

%------------------------------------------------------------------------------%
% Initialize
%--------------------------------------------------------------------------
nloop = 0;
loop = 1;
status=0;
goldensearch_deltmaxinit= option.goldensearch_deltmax;




%-----------------------------------------
% Initializing SVM
%------------------------------------------
SumSigma=sum(Sigma);
[radius, optimalBeta,xsup,Alpsup,w0,pos,alpha,margin,obj] = radiusMarginBoundClass([],yapp,lambdareg,K,Sigma,verbosesvm,alphainit);
%compute gradient
[grad] = radiusMarginBoundGrad(K,Sigma,pos,Alpsup,radius,optimalBeta);



Sigmaold  = Sigma ;
history.obj=[];
history.sigma=[];
history.KKTconstraint=[1];
history.dualitygap=[1];

history.iteration=0;
%------------------------------------------------------------------------------%
% Update Main loop
%------------------------------------------------------------------------------%

while loop & nloopmax >0 ;

    history.iteration = nloop;
    
    nloop = nloop+1;
    history.sigma= [history.sigma;Sigma];
    history.obj=[history.obj obj];

    %-----------------------------------------
    % Update weigths Sigma
    %-----------------------------------------
    t = cputime ;
    [radius,optimalBeta, margin,Sigma,Alpsup,w0,pos,obj] = radiusMarginBoundUpdate(K,Sigma,pos,Alpsup,w0,yapp,grad,obj,radius,optimalBeta, margin, option) ;
    
    %-----------------------------------------
    % Thresholding
    %-----------------------------------------

%     %-------------------------------
%     % Numerical cleaning
%     %-------------------------------
    Sigma(find(abs(Sigma<option.numericalprecision)))=0;
    Sigma=Sigma/sum(Sigma);
    %-----------------------------------------------------------
    % Enhance accuracy of line search if necessary
    %-----------------------------------------------------------
    if max(abs(Sigma-Sigmaold))<option.numericalprecision & option.goldensearch_deltmax > optiongoldensearchmax
        option.goldensearch_deltmax=option.goldensearch_deltmax/10;
    elseif option.goldensearch_deltmax~=goldensearch_deltmaxinit
        option.goldensearch_deltmax*10;
    end;

    %----------------------------------------------------
    % process approximate KKT conditions
    %--------------------------------------------------
    [grad] = radiusMarginBoundGrad(K,Sigma,pos,Alpsup,radius,optimalBeta);

    indzero=find(abs(Sigma<=option.numericalprecision));

    %----------------------------------------------------
    % process duality gap
    %----------------------------------------------------
    normek=-grad ; % 0.5*Alpsup'*K(pos,pos,i)*Alpsup;  

    %------------------------------------------
    %  verbosity
    %------------------------------------------
    if verbose
        if nloop == 1 || rem(nloop,10)==0
            fprintf('--------------------------------------------------\n');
            fprintf('Iter | Obj.    | DiffBetas\n');
            fprintf('--------------------------------------------------\n');
        end;
        fprintf('%d   | %8.4f | %6.4f \n',[nloop obj   max(abs(Sigma(2:end)-Sigmaold(2:end)))]);
    end


    %----------------------------------------------------
    % check variation of Sigma conditions
    %----------------------------------------------------
    if  option.stopvariation==1 & option.stopKKT== 0 & max(abs(Sigma(2:end) - Sigmaold(2:end)))<option.seuildiffsigma;
        loop = 0;
        fprintf(1,'variation convergence criteria reached \n');
        history.sigma= [history.sigma;Sigma];
        history.obj=[history.obj obj];
    end;

    %-----------------------------------------------------
    % check nbiteration conditions
    %----------------------------------------------------
    if nloop>=nloopmax ,
        loop = 0;
        history.sigma= [history.sigma;Sigma];
        history.obj=[history.obj obj];
        status=2;
        fprintf(1,'maximum number of iterations reached\n')
    end;
    if nloop < option.miniter & loop==0
        loop=1;
    end;

    %-----------------------------------------------------
    % Updating Variables
    %----------------------------------------------------

    Sigmaold  = Sigma ;
end;
