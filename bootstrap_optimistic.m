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
    [f, itn] = newton_from_zero(alpha, v, R, tol, maxit);
    if sum(f) > 1 - sqrt(eps)
        warning 'is your problem supercritical?'
    end
    total_iterations = total_iterations + itn;
    beta = 1 - sum(f);
    y = x-f;
    y = y / sum(y) * beta;
    for it = 1:maxit-total_iterations
        matrix = alpha*R*(kron(f+y,eye(n)) + kron(eye(n),f));
        [oldy, y] = deal(y, perronvector(matrix,'eig'));
        y = y / sum(y) * beta;
        x = y + f; res = alpha*R*kron(x, x) + (1-alpha)*v - x;
        if norm(res, 1) < tol
            break
        end
    end
    total_iterations = total_iterations + it;
    if total_iterations >= maxit
        warning('Maximum numer of iterations reached');
    end
end
it = total_iterations;