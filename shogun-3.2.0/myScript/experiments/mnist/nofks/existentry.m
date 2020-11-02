function b = existentry(M,r,c)

if (size(M,1)<r) 
    b = 0;
    return;
end;

if (size(M,2)<c)
    b = 0;
    return;
end;

if(M(r,c)~=0)
    b = 1;
    return;
end;

b = 0;
return;



