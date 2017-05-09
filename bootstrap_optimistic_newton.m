function [x, it] = bootstrap_optimistic_newton(alpha, v, R, tol, maxit)
% Successive Newtons incrasing the value of alpha

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

if alpha > 0.6
    alphas = linspace(0.6, alpha, 20);
else
    alphas = alpha;
end

n = length(v);
x = v;

total_iterations = 0;

for alpha = alphas
    [x, it] = optimistic_newton(alpha, v, R, tol, maxit-total_iterations, x);
    total_iterations = total_iterations + it;
end
it = total_iterations;