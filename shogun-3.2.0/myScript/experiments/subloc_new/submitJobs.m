
%datasetName = 'plant';
datasetName = 'nonpl';
%datasetName = 'psortNeg';
%datasetName = 'psortPos';

%optName = 'DIRECT';
%optName = 'NEWTON1';
optName = 'CPLEX';

epsSetting = 'b';
splits = 1:30;

cmdCplexPath = 'export ILOG_LICENSE_FILE=/fml/ag-raetsch/share/software/ilog/licenses/access-587380.ilm';
%cmdCplexPath = 'export ILOG_LICENSE_FILE=/fml/ag-raetsch/share/software/ilog/licenses/access-703040.ilm';
%cmdCplexPath = 'export ILOG_LICENSE_FILE=/fml/ag-raetsch/share/software/ilog/licenses/access-820980.ilm';



% === prep

cmdCd = 'cd ~/work/MKL-q/MKL/lp_mkl/subloc/';
cmdMatlab = 'matlab -nojvm -r "runit( ''%s'', ''%s'', ''%s'', %d ); quit"';

switch( optName )
 case 'CPLEX'
  %cmdShogunPath = 'export LD_LIBRARY_PATH=~/shogun/branches/mkl/src/libshogun:~/shogun/branches/mkl/src/libshogunui';
  cmdShogunPath = 'export LD_LIBRARY_PATH=~/shogun/branches/mkl/src/MYLIB';
  optNameShort = 'C';
 case 'DIRECT'
  cmdShogunPath = 'export LD_LIBRARY_PATH=~/shogun/branches/mkl/src/MYLIB';
  optNameShort = 'D';
 case { 'NEWTON1' }
  cmdShogunPath = 'export LD_LIBRARY_PATH=~/shogun/branches/mkl/src/MYLIB';
  optNameShort = [ 'N' optName(end) ];
 otherwise
  error( '???' );
end;

switch( datasetName )
 case 'plant'
  datasetNameShort = 'pl';
  reqMem = 2;
 case 'nonpl'
  datasetNameShort = 'np';
  reqMem = 7;
 case 'psortPos'
  datasetNameShort = 'P+';
  reqMem = 1;
 case 'psortNeg'
  datasetNameShort = 'P-';
  reqMem = 2.5;
 otherwise
  error( '???' );
end;

cplexLicence = cmdCplexPath( length('export ILOG_LICENSE_FILE=!'):end );
ASSERT( cplexLicence((end-3):end) == '.ilm' );
ASSERT( cplexLicence(1:52) == '/fml/ag-raetsch/share/software/ilog/licenses/access-' );
cplexLicence = cplexLicence( 53:58 );
switch( optName )
 case 'CPLEX'
  resDir = [ 'CPLEX-' cplexLicence '_' epsSetting ];
 otherwise
  resDir = [ optName '_' epsSetting ];
end;
resDir = [ 'res/' datasetName '/' resDir '/' ];



% === submit

for( numPerm = splits )
  jobName = sprintf( '%s%02d%s%s', datasetNameShort, numPerm, optNameShort, epsSetting );
  
  % --- computations already done?
  resFileName = sprintf( '%s%s,split_%02d.mat', resDir, datasetName, numPerm );
  if( exist(resFileName,'file') )
    savedVars = who( '-FILE', resFileName );
    if( length(savedVars) < 3 )
      fprintf( '# %s: probably corrupt\n', jobName, unfinished );
    end;
  end;
  if( exist(resFileName,'file') )
    load( resFileName, 'RES' );
    R = RES( :, :, 1 );
    r = R(:);
    unfinished = sum(r<0) + sum(isnan(r));
    if( unfinished == 0 )
      fprintf( '# %s: finished\n', jobName );
      continue;
    end;
    %fprintf( '# %s: %d unfinished\n', jobName, unfinished );
  end;
  
  % --- job script
  jobFileName = sprintf( 'jobs/%s.sh', jobName );
  fid = fopen( jobFileName, 'w' );
  if( ~strcmp(optName,'DIRECT') )
    fprintf( fid, '%s\n', cmdCplexPath );
  end;
  fprintf( fid, '%s\n', cmdShogunPath );
  fprintf( fid, '%s\n', cmdCd );
  fprintf( fid, [ cmdMatlab '\n' ], datasetName, optName, epsSetting, numPerm );
  fclose( fid );
  
  % --- qsub
  qsub = sprintf( 'qsub  -l h_vmem=%.1fG  -N %s  -o $HOME/tmp/log/  -e $HOME/tmp/log/  jobs/%s.sh', reqMem, jobName, jobName );
  fprintf( [ qsub '\n' ] );
end;



