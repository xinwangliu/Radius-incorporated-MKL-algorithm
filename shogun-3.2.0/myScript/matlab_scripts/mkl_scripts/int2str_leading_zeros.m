function str = int2str_leading_zeros( d, i)
% d  integer
% i   leading zero offset

istr = int2str(d);
str = [char(48*ones(1,i-length(istr))), istr];
