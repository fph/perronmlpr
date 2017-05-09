function [x, it] = bootstrap_optimistic_newton(alpha, v, R, tol, maxit)
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

total_iterations = 0;
x = v;
for i = 1:length(alphas)
    alpha = alphas(i);
    if i == 1
        x_guess = v;
    else
        % updating x with a first-order estimate
        fx = alpha*R*kron(eye(n),x) + alpha*R*kron(x,eye(n)) - eye(n);
        falpha = R*kron(x,x) - v;
        x_guess = x + (alpha - alphas(i-1)) * (-fx \ falpha);
    end
    x_old = x;
    [x, it] = optimistic_newton(alpha, v, R, tol, maxit-total_iterations, x_old);
    if i>1
        disp('------');
        second_order_deriv = norm(x_guess - x) * 2 / norm(alpha - alphas(i-1))^2
        norm(x_guess - x)
        norm(x_old - x)
    end
    total_iterations = total_iterations + it;
end
it = total_iterations;