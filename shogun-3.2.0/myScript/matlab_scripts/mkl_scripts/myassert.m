
function assert( assertions, message );

if( ~exist('assertions','var') )
  error( 'assert: not enough input arguments.' );
end;
if( all( assertions(:) ) )
   % if assertions are true, do nothing
   return;
end
if( ~exist('message','var') )
   message = 'Assertion failure';
end

error( message );

