function [x, it, itn, xs] = optimistic_fixed_point(alpha, v, R, tol, maxit);
% new fixed point version of the "optimistic iteration"
% starts from x0 = 0

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end
n = length(v);
scaledR = alpha*R;

xs = [];

[f, itn] = newton_from_zero(alpha, v, R);
beta = 1 - sum(f);
y = beta * ones(n,1)/n;
for it = 1:maxit-itn
    constant = scaledR*(kron(f,y) + kron(y,f));
    matrix = scaledR*kron(y,eye(n));
    [oldy, y] = deal(y, (eye(n)-matrix) \ constant);
    y = y / sum(y) * beta;
    x = y + f; res = alpha*R*kron(x, x) + (1-alpha)*v - x;
    if norm(res, 1) < tol
        break
    end
    if nargout > 3
        xs = [xs f + y];
    end 
end
x = y + f;
if it == maxit
    warning('Maximum numer of iterations reached');
end