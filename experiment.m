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

xs = [];

total_iterations = 0;

R = 1/2*(R + R*commutation(n,n));

for alpha = alphas
    [f, itn] = newton_from_zero(alpha, v, R, tol, maxit);
    if sum(f) > 1 - sqrt(eps)
        warning 'is your problem supercritical?'
    end
    total_iterations = total_iterations + itn;
    beta = 2- 1/ alpha;
    y = x-f;
    y = y / sum(y) * beta;
    for it = 1:maxit-total_iterations
        matrix = alpha*R*(kron(f+y,eye(n)) + kron(eye(n),f));
        [u, lambda] = perronvector(matrix, 'eig'); %lambda = 1 should hold
        u = u / sum(u) * beta;
        e = ones(n,1);
        Jac = (eye(n)-matrix+u*e'/(e'*u)) \ (alpha*R*kron(eye(n),u)) - alpha*u*e';
        [oldy, y] = deal(y, y - (eye(n) - Jac) \ (y - u));
        y(y<0) = 0; 
        y = y / sum(y) * beta;
        x = y + f; res = alpha*R*kron(x, x) + (1-alpha)*v - x;
        if norm(res, 1) < tol
            break
        end
        xs = [xs x];
    end
    total_iterations = total_iterations + it;
    if total_iterations >= maxit
        warning('Maximum numer of iterations reached');
    end
end
it = total_iterations;
