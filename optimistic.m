function [x, it, itn, xs, res, rho_jac] = optimistic(alpha, v, R, tol, maxit);
% new Perron vector method
% starts from x0 = 0

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end
n = length(v);
scaledR = alpha*R;
scaledv = (1-alpha)*v;

xs = [];
res = [];

[f, itn] = newton_from_zero(alpha, v, R, tol, maxit);
if sum(f) > 1 - sqrt(eps)
    warning 'is your problem supercritical?'
end
beta = 1 - sum(f);
y = beta * ones(n,1)/n;
for it = 1:maxit-itn
    matrix = scaledR*(kron(f+y,eye(n)) + kron(eye(n),f));
    [oldy, y] = deal(y, perronvector(matrix,'eig'));
    y = y / sum(y) * beta;
    x = y + f; res = alpha*R*kron(x, x) + (1-alpha)*v - x;
    if norm(res, 1) < tol
        break
    end
    if nargout > 3
        xs = [xs f + y];
        x = y + f;
        res = [res norm(x - scaledR*kron(x,x) - scaledv)];
    end
end
x = y + f;
it = itn + it;
if it == maxit
    warning('Maximum numer of iterations reached');
end
if nargout > 3
    e=ones(n,1);
    jac=(eye(n)- y*e'/(e'*y)) * pinv(matrix - eye(n)) * (y*e'/(e'*y) - eye(n)) * scaledR * kron(eye(n),y);
    rho_jac=max(abs(eig(jac)));
end