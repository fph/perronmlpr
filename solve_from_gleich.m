function [xn, it] = shifted(alpha, v, R, tol, maxit, x0)
% emulates solve() in Gleich's code
% starts from x0 a stochastic vector
% Contains Gleich's code without the OO part

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

if alpha< 1/2
    [xn, it] = shifted(alpha, v, R, tol, maxit);
else
    [xn, it] = innout(alpha, v, R, tol, maxit);
end