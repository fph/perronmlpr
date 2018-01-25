function [x, it] = bootstrap_wrong_derivative(inner_solver, target_alpha, v, R, tol, maxit, relative_speed)
% Calls an inner solver iteratively over increasing values of alpha
% Predicts the new x using a first-order Taylor expansion with a "modified"
% (i.e., wrong) first derivative

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end
if not(exist('relative_speed', 'var')) || isempty(relative_speed)
    relative_speed = 0.01;
end

n = length(v);

total_iterations = 0;
old_x = nan;
old_alpha = nan;

alpha = 0.6;
while true
    if any(isnan(old_x))
        x_guess = v;
    else
        % updating x with a first-order estimate, with a "wrong" choice of
        % alpha (should have used old_alpha; instead we use alpha in the
        % new point).
        fx = alpha*R*kron(eye(n),old_x) + alpha*R*kron(old_x,eye(n)) - eye(n);
        falpha = R*kron(old_x,old_x) - v;
        x_guess = old_x + (alpha - old_alpha) * (-fx \ falpha);
    end
    [x, it] = inner_solver(alpha, v, R, tol, maxit-total_iterations, x_guess);
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
        second_derivative_guess = norm(x_guess - x) * 2 / norm(alpha - old_alpha)^2;
        step_size = sqrt(2*relative_speed / second_derivative_guess);
        new_alpha = alpha + step_size;
    end
    if new_alpha > target_alpha
        new_alpha = target_alpha;
    end
%    alpha, new_alpha
%    [x_guess-x]
%      xguessdiff = norm(x_guess-x) %TODO: debug print
%      new_alpha %TODO: debug print
    [alpha, old_alpha] = deal(new_alpha, alpha);
    old_x = x;
    if total_iterations >= maxit
        break
    end
end
it = total_iterations;