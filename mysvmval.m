function [y] = mysvmval(w,b,pos,Ktsttrn)
y = Ktsttrn(:,pos)*w + b;