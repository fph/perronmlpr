function [x, it] = newton_from_zero(alpha, v, R, tol, maxit);
% Newton iteration in "Multilinear Pagerank"
% starts from x0 = 0

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

n = length(v);
scaledv = (1-alpha)*v;
scaledR = alpha*R;
x = zeros(n,1);
for it = 1:maxit
    res = scaledR*kron(x,x) + scaledv - x;
    Jac = eye(n) - scaledR * (kron(x,eye(n)) + kron(eye(n),x));
    [oldx, x] = deal(x, x + Jac\res);
    if norm(oldx-x) < tol
        break
    end
end
if it == maxit
    warning('Maximum numer of iterations reached');
end