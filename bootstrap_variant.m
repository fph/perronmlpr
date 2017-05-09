function [x, it] = bootstrap(alpha, v, R, tol, maxit)
% Successive Newtons incrasing the value of alpha

target_alpha = alpha;
n = length(v);

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

assert(target_alpha >= 0.7); %TODO: for now we focus on this case
alpha = 0.7;
[x, it] = fixed_point(alpha, v, R);
total_iterations = it;

while alpha < target_alpha
    
    % selects epsilon as large as possible such that
    % the fixed point eqn is a contraction, i.e., 4*norm(RR)*norm(vv) < 1
    epsilon = 0.0000001;
    M = eye(n) - (alpha+epsilon)*(R*kron(eye(n),x)+R*kron(x,eye(n)));
    vv = epsilon*(M\(R*kron(x,x) - v));
    RR = (alpha + epsilon) * (M\R);
    assert(4*norm(RR)*norm(vv)<1);
    
    while 4*norm(RR)*norm(vv)<1
        last_good_epsilon = epsilon;
        epsilon = epsilon * 1.1;
        M = eye(n) - (alpha+epsilon)*(R*kron(eye(n),x)+R*kron(x,eye(n)));
        vv = epsilon*(M\(R*kron(x,x) - v));
        RR = (alpha + epsilon) * (M\R);
    end
%    4*norm(RR)*norm(vv)
    epsilon = last_good_epsilon;
    
    if alpha + epsilon >= target_alpha
        epsilon = target_alpha - alpha;
    end
    
    alpha = alpha + epsilon;
    [x, it] = newton(alpha, v, R, tol, maxit-total_iterations, x);
    total_iterations = total_iterations + it;
end