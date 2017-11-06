function [x, it] = bootstrap_optimistic_newton(target_alpha, v, R, tol, maxit)
% Successive Newtons incrasing the value of alpha

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

n = length(v);
relative_speed = 0.01;

total_iterations = 0;
old_x = nan;
old_alpha = nan;

alpha = 0.6;
while true
    if any(isnan(old_x))
        x_guess_firstorder = nan;
        x_guess_secondorder = nan;
        x_guess_secondorder_bea = nan;
        x_guess = v;
    else
        % updating x with a second-order estimate
        % partial derivatives of alpha*R*kron(x,x)+(1-alpha)*v-x
        partialx = old_alpha*R*kron(eye(n),old_x) + old_alpha*R*kron(old_x,eye(n)) - eye(n);
        partialalpha = R*kron(old_x,old_x) - v;
        % from the implicit function theorem
        xprime = -partialx \ partialalpha;
        % now we differentiate (in alpha, total derivative) this expression for xprime piece by piece.
        partialalpha_prime = R*kron(old_x, xprime) + R*kron(xprime, old_x);
        partialx_prime = R*kron(eye(n),old_x) + R*kron(old_x,eye(n)) + old_alpha*R*kron(eye(n), xprime) + old_alpha*R*kron(xprime, eye(n));
        % derivative of inv(partialx) =
        % -inv(partialx)*partialx_prime*inv(partialx)
        xsecond = -partialx \ (partialx_prime * xprime) - partialx \ partialalpha_prime;
        x_guess_firstorder = old_x + (alpha - old_alpha) * xprime;
        x_guess_secondorder = old_x + (alpha - old_alpha) * xprime + 1/2*(alpha - old_alpha)^2 * xsecond;
        x_guess = x_guess_secondorder;
    end
    [x, it] = optimistic_newton(alpha, v, R, tol, maxit-total_iterations, x_guess);
    total_iterations = total_iterations + it;
    if alpha >= target_alpha
        break
    end
    % construct new alpha
    if any(isnan(old_x))
        % at the first step, we have no "second derivative" information
        % available
        new_alpha = alpha + 0.01;
    else
        third_derivative_guess = norm(x_guess - x) * 6 / norm(alpha - old_alpha)^3;
        step_size = (6*relative_speed / third_derivative_guess)^(1/3);
        new_alpha = alpha + step_size;
    end
    if new_alpha > target_alpha
        new_alpha = target_alpha;
    end
%     xguessdiff = norm(x_guess-x) %TODO: debug print
%     new_alpha %TODO: debug print
    [alpha, old_alpha] = deal(new_alpha, alpha);
    old_x = x;
    if total_iterations >= maxit
        break
    end
end
it = total_iterations;