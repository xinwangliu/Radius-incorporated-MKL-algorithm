
function [X,y] = mnist_load(nofex, rnd_seed)


% load data
load ../data/mnist;

% seed random number generator
rand('seed',rnd_seed);

% shuffle
idx = randperm(size(Xtr,1));

% draw
X = Xtr(idx(1:nofex),:);
y = ytr(idx(1:nofex));

% transform
X = X./255;
y = mod(y,2);
y(y==0) = -1;

clear Xtr;
clear ytr;
clear Xte;
clear yte;

return;

