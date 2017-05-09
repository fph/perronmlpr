function [x it] = tpr_method(method, alpha, v, R, tol, maxit)
% Runs one of the methods in Gleich's tensorpr3 class
%
% result = tpr_method(alpha, v, R, method, tol, maxit)

addpath('./mlpagerank-master')

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

tpr = tensorpr3(R, alpha, v);
[x hist] = tpr.(method)('maxiter',maxit, 'tol', tol);
it = length(hist);