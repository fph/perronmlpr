function [xn, it] = innout(alpha, v, R, tol, maxit, x0)
% Inner-outer iteration
% starts from x0 a stochastic vector
% Contains Gleich's code without the OO part

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

n = length(v);

% if not(exist('x0','var')) || isempty(x0)
%     x0 = (1-alpha)*v;
% end


Rt = alpha*R + (1-alpha)*v*ones(1, n^2);
at = alpha / 2;
x = v;

for it = 1:maxit
    xn = solve_from_gleich(at, x, Rt, max(tol/10,eps(1))); %TODO: count inner iterations
    xn = xn/sum(xn);
    
    res = alpha*R*kron(xn, xn) + (1-alpha)*v - xn;
    %curdiff = norm(x - xn,1);

    % check for termination
    if norm(res, 1) <= tol %|| curdiff <= tol/10 %TODO: remove this
        %norm(res, 1)
        break;
    end
    x = xn;
end
