function [ betas, classy, ResMKL ] = train_mklSVM( Kxs, y, varargin )

% train_mklSVM - trains a MKL - Suport Vector Machine (SVM) 
%
% Synopsis:
%   [betas, classy]= train_mklSVM( Kxs, y, properties)
%

%==============================================================================
% INIT
%==============================================================================

paraSVM= propertylist2struct(varargin{:});
paraSVM= set_defaults(paraSVM, ...
                         'implementation', 'LIBSVM', ...
                         'C',1, ...
                         'verbosity', 0, ...
                         'epsilon', 10^(-3), ...
                         'use_mosek', 0 ,...
                         'addRidge', 0 );
PARA.verbosity = paraSVM.verbosity;
                         
% === init
nofKernels = size( Kxs, 3 );
mAll = size( Kxs,1 );

% === check
% --- dims
myassert( size(Kxs,1) == mAll );
myassert( size(Kxs,2) == mAll );
% --- content
myassert( isfinite(Kxs) );

% === mklParams
if PARA.verbosity==0
  warning off;
end
% --- model
PARA.maxBeta = 1;
% --- optimization
PARA.epsRidge = 1e-10;
PARA.epsBeta = 1e-12;
% --- convergence criteria
PARA.deltaObj = paraSVM.epsilon;
PARA.deltaBetas = 0;
PARA.maxinfeasibility = 1e-6;

%==============================================================================
% Prepare SILP
%==============================================================================

% === auxilliary
theta = 0;
betas = ones( nofKernels, 1 ) / nofKernels;
oldTheta = nan;
oldBetas = nan * betas;
% --- kernels
if paraSVM.addRidge
  for( k = 1:nofKernels )
    Kxs(:,:,k) = Kxs(:,:,k) + PARA.epsRidge * mean(diag(Kxs(:,:,k))) * speye(mAll);
    if PARA.verbosity
      fprintf( 'kernel: %d\n', k);
    end
  end
end
if PARA.verbosity
  fprintf( '\n' );
end

% === betas (LP) [grows during iterations]
SILP = [];
SILP.verbosity = 0;
% --- objective
SILP.f = zeros( nofKernels+1, 1 );
SILP.f(end) = 1;  % the last variable is theta
% --- equality constraints
SILP.A = ones( 1, nofKernels+1 );
SILP.A(end) = 0;
SILP.b = 1;
% --- bounds
SILP.LB = zeros( nofKernels+1, 1 );
SILP.LB(end) = -inf;
SILP.UB = PARA.maxBeta * ones( nofKernels+1, 1 );
SILP.UB(end) = +inf;


%==============================================================================
% Solve SILP
%==============================================================================

% === iterate optimizing alphas / betas
iter = 0;
ResMKL = [];
ResMKL.thetas = [];
ResMKL.betas = {};
ResMKL.SVM = {};
while( 1 )
  iter = iter + 1;
  ResMKL.betas{iter} = betas;
  
  % === base learner
  K = reshape( reshape( Kxs, mAll*mAll, nofKernels) * betas  ,  mAll,mAll, 1);
  ResSVM = train_sgSVM( 'y', y, 'Ktr', K,'C', paraSVM.C, 'implementation', paraSVM.implementation);
  alpha = zeros(mAll,1);
  alpha(ResSVM.alphas(:,2)+1) = ResSVM.alphas(:,1);
  ResMKL.SVM{iter} = ResSVM;

  % === convergence?
  theta = ResSVM.objVal;
  %ResMKL.thetas(end+1) = theta;
  %dTheta = abs( theta - oldTheta );
  % Ulf's Abbruchbedingung leicht geändert!
  dTheta = abs(1-  theta/ oldTheta );
  dBetas = norm( betas - oldBetas, 'fro' );
  % --- output
  if( PARA.verbosity >= 1 )
    fprintf( 'SILP iter %2d: theta=%.4f, d_theta=%.4f, d_B=%.4f\n', iter, theta, dTheta, dBetas );
    %disp(betas);
  end;
  % --- convergence check
  if( nofKernels == 1 )
    break;
  end;
  if( iter > 1 )
    if( dTheta < PARA.deltaObj )
      break;
    end;
    if( dBetas < PARA.deltaBetas )
      break;
    end;
  end;
  % === betas (LP)
  linC = zeros( nofKernels, 1 );
  for( k = 1:nofKernels )
    linC(k) = sum(abs(alpha)) - 0.5* alpha'*Kxs(:,:,k)*alpha;
     %Achtung in \alpha sind schon die labels drin! alpha_i --> alpha_i*y_i
  end;
  myassert( size(linC) == [nofKernels,1] );
  SILP.A( end+1, : ) = [ linC' -1 ];
  SILP.b( end+1, 1 ) = 0;
  
  oldTheta = theta;
  oldBetas = betas;
  x0 = [ betas(:) ; 0 ];
  t = SILP.A*x0 - SILP.b;
  x0(end) = max(t);
  t = SILP.A*x0 - SILP.b;
  if( max(t) >= 1e-12 )&&PARA.verbosity
    fprintf( 'WARNING: numerical infeasibility = %e\n', max(t) );
  end;
  myassert( t <= PARA.maxinfeasibility );
   
  options = optimset( 'Display', 'off', 'LargeScale', 'off' );
  % --- optimize SILP
  if paraSVM.use_mosek
    clear prob;
    prob.c = SILP.f;
    prob.a = sparse(SILP.A);
    prob.buc = SILP.b;
    prob.blc = [SILP.b(1); -inf * ones(size(SILP.b,1)-1,1)];
    prob.blx = SILP.LB;
    prob.bux = SILP.UB;
    %prob.sol = x0;
    if PARA.verbosity==2
      [flag, res] = mosekopt('minimize echo(1)', prob);
    else
      [flag, res] = mosekopt('minimize echo(0)', prob);
    end
    if exist('res') && isfield(res,'sol') && isfield(res.sol,'itr') && isfield(res.sol.itr,'prosta') && equal(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
      x0 = res.sol.itr.xx;
      fval = SILP.f'*x0;
      flag = 0;
    else
      flag = -1;
      x0 = flipud(x0);
      % quasi neustarten, mit zufälligem Element
    end
  else
    [ x0, fval, flag, output, lambda ] = linprog( SILP.f, SILP.A(2:end,:), SILP.b(2:end), SILP.A(1,:), SILP.b(1), SILP.LB, SILP.UB, x0, options );
  end
      
   % assert( flag >= 0 );
  ResSILP.val = fval;
  % ResSILP.lambda = lambda;

  oldtheta = x0(end);
  ResMKL.thetas(end+1) = theta;
  betas = reshape( x0(1:(end-1)), nofKernels, 1 );
  betas( betas < PARA.epsBeta ) = 0;
  if( PARA.verbosity >= 2 )
    betaStr = sprintf( '%6.4f ', betas );
    fprintf( '  betas = %s\n', betaStr );
  end;
  
end;


%==============================================================================
% DONE
%==============================================================================

% === output
if( PARA.verbosity == 1 )
  betaStr = sprintf( '%6.4f ', betas );
  fprintf( '  betas = %s\n', betaStr );
end;

classy = ResSVM;