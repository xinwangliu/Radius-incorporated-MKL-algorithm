
function [ K ] = rbf( X, sigma )
N = size(X,1);
K = (repmat(sum((X.^2)', 1), [N 1])' + repmat(sum((X.^2)', 1), [N 1]) - 2*X*(X'));
K = exp(-K./sigma);
return;
