
function [ opt ] = setDefault( opt, param, value );
  
if( ~ isfield(opt,param) )
  opt = setfield( opt, param, value );
end;


