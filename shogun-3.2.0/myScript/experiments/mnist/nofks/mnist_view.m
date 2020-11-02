
function mnist_view(vec,label)

% show label
fprintf('label = %d\n',label);

% print image
for i=1:28
    for j=1:28
        fprintf('%4d ',vec((i-1)*28+j));
    end;
    fprintf('\n');
end;

        