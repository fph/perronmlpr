function [x, it, itn] = newton_on_y(alpha, v, R, tol, maxit)
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

[f, itn] = newton_from_zero(alpha, v, R);
beta = 1 - sum(f);
y = beta * ones(n,1)/n;

for it = 1:maxit
    res = scaledR*(kron(f+y,y)+kron(y,f)) - y;
    Jac = eye(n) - scaledR * (kron(y + f,eye(n)) + kron(eye(n),y + f));
    [oldy, y] = deal(y, y + Jac\res);
    y = y / sum(y) * beta;
    if norm(oldy-y) < tol
        break
    end
end
if it == maxit
    warning('Maximum numer of iterations reached');
end
x = y + f;