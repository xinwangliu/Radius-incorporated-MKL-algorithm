
function [ RES, INFO ] = playNonpl( numPerm );



% === Init

% dbstop if error;

% --- start shogun
addpath( '~/shogun/branches/mkl/src/matlab/' );
version = sg( 'get_version' );

% --- select available cplex license server
addpath( '/fml/ag-raetsch/home/zien/svn/tools/cplex_matlab-7.0/cplex91/' );
cplex_license;



% === Param

param0 = [];
param0.datasetName = 'nonpl';
param0.numPerm = numPerm;

% --- setting
param0.fracTst = 0.2;
param0.nofFoldsVal = 3;
param0.foldVal = 0;

% --- model
param0.mcMode = '1vsRest';
%param0.mklSolver = 'INTERNAL';
param0.mklSolver = 'CPLEX';
param0.logCs = -5:2:+5;
param0.mklNorms = [ 1 8/7 4/3 2 4 8 ];

% --- technical
param0.cacheSize = 10;
param0.svmEps = 1e-3;
param0.mklEps = 1e-4;



% === Experiments

[ RES, INFO ] = trainAndTestSeveral( param0 );


