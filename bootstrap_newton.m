function [x, it] = bootstrap(alpha, v, R, tol, maxit)
% Successive Newtons incrasing the value of alpha

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

if alpha > 0.7
    alphas = linspace(0.7, alpha, 20);
else
    alphas = alpha;
end

n = length(v);
x = v;
total_iterations = 0;
for alpha = alphas
    for it = 1:maxit-total_iterations
        scaledv = (1-alpha)*v;
        scaledR = alpha*R;
        res = scaledR*kron(x,x) + scaledv - x;
        Jac = eye(n) - scaledR * (kron(x,eye(n)) + kron(eye(n),x));
        [oldx, x] = deal(x, x + Jac\res);
        res = alpha*R*kron(x, x) + (1-alpha)*v - x;
        if norm(res, 1) < tol
            break
        end
        x = x / sum(x);
    end
    total_iterations = total_iterations + it;
    if total_iterations >= maxit
        warning('Maximum numer of iterations reached');
    end
end
it = total_iterations;