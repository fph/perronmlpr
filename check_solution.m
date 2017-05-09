function result = check_solution(alpha, v, R, x, tol)
% checks the validity of a solution to the MLPR equation and returns an
% error code
%
% check_solution(alpha, v, R, x)
% 
% 0 = all good
% 1 = residual is not zero numerically
% 2 = solution is not real (numerically)
% 3 = solution is not positive (numerically)
% 4 = solution has not mass 1

n = length(v);

residual = alpha * R * kron(x, x) + (1-alpha) * v - x;

if norm(residual, 1) / norm(x, 1) > tol
    result = 1; return
end

if imag(x) / norm(x, 1) > tol
    result = 2; return
end

if min(x) / norm(x, 1) < - tol
    result = 3; return
end

if abs(sum(x) - 1) > tol
    result = 4; return
end

result = 0;