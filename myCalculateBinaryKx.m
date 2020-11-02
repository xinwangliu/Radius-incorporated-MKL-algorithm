function [Kx] = myCalculateBinaryKx(w,pos,num)

Kx = zeros(num);
Kx(pos,pos) = (w/norm(w))*(w'/norm(w));