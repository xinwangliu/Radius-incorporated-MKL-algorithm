
function [ RES, INFO ] = playPlant( optName, numPerm );



% === Init

% dbstop if error;

% --- const
epsSetting = 'a';

% --- start shogun
switch( optName )
 case 'CPLEX'
  mklSolver = 'CPLEX';
  %addpath( '~/shogun/branches/mkl/src/matlab/' );
  addpath( [ '~/shogun/branches/mkl/src/MYSG/NEWTON1/' ] );
 case { 'DIRECT', 'NEWTON1' }
  mklSolver = 'INTERNAL';
  addpath( [ '~/shogun/branches/mkl/src/MYSG/' optName '/' ] );
  %%%addpath( '~/shogun/branches/mkl/src/matlab/' );
 otherwise
  error( '???' );
end;
sgVersion = sg( 'get_version' )
fprintf( '\n' );

% --- select available cplex license server
if( strcmp(mklSolver,'CPLEX') )
  if( 0 )
    addpath( '/fml/ag-raetsch/home/zien/svn/tools/cplex_matlab-7.0/cplex91/' );
    cplex_license;
  end;
  cplex = getenv( 'ILOG_LICENSE_FILE' );
  ASSERT( cplex((end-3):end) == '.ilm' );
  ASSERT( cplex(1:52) == '/fml/ag-raetsch/share/software/ilog/licenses/access-' );
  cplex = cplex( 53:58 );
  fprintf( 'CPLEX license: %s\n', cplex );
end;

% --- res dir
switch( mklSolver )
 case 'CPLEX'
  resDir = [ 'CPLEX-' cplex '_' epsSetting ];
 case 'INTERNAL'
  resDir = [ optName '_' epsSetting ];
 otherwise
  error( '???' );
end;
success = mkdir( 'res/', resDir );
ASSERT( success );
resDir = [ 'res/' resDir '/' ];
fprintf( 'results -> %s\n', resDir );
fprintf( '\n' );



% === Param

param0 = [];
param0.datasetName = 'plant';
param0.numPerm = numPerm;
param0.resDir = resDir;

% --- setting
param0.fracTst = 0.2;
param0.nofFoldsVal = 3;
param0.foldVal = 0;

% --- model
param0.optName = optName;
param0.mklSolver = mklSolver;
param0.mcMode = '1vsRest';
param0.logCs = [ -5 -3 -1 0 +1 +2 +3 +5 +7 ];
%param0.logCs = [ -5 -3 -1 +1 +2 +3 +5 +7 ];
param0.mklNorms = [ 1 32/31 16/15 8/7 4/3 2 4 8 16 +inf ];
%param0.mklNorms = [ 1 32/31 16/15  8 16 +inf ];
%param0.mklNorms = [ 4/3 2 4 ];

% --- technical
param0.cacheSize = 10;
param0.epsSetting = epsSetting;
switch( param0.epsSetting )
 case 'a'
  param0.svmEps = 1e-3;
  param0.mklEps = 1e-5;
 case 'b'
  param0.svmEps = 1e-2;
  param0.mklEps = 1e-4;
 otherwise
  error( '???' );
end;



% === Experiments

[ RES, INFO ] = trainAndTestSeveral( param0 );


