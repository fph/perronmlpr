function [x, it] = fixed_point(alpha, v, R, tol, maxit, x0);
% fixed point iteration in "Multilinear Pagerank"
% starts from x0 = v unless specified

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end
if not(exist('x0','var')) || isempty(x0)
    x0 = v;
end

x = x0;
scaledv = (1-alpha)*v;
scaledR = alpha*R;
for it = 1:maxit
    [oldx, x] = deal(x, scaledR*kron(x,x) + scaledv);
    res = alpha*R*kron(x, x) + (1-alpha)*v - x;
    if norm(res, 1) < tol
        break
    end
    x = x / sum(x);
end
if it == maxit
    warning('Maximum numer of iterations reached');
end