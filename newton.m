function [xn, it] = newton(alpha, v, R, tol, maxit, x0)
% Newton iteration in "Multilinear Pagerank"
% starts from x0 a stochastic vector
% Contains Gleich's code without the OO part

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

n = length(v);
scaledv = (1-alpha)*v;
scaledR = alpha*R;

if not(exist('x0','var')) || isempty(x0)
    x0 = (1-alpha)*v;
end

xcur = x0;
for it = 1:maxit
    I = eye(n);
    A = alpha*R*(kron(xcur, I) + kron(I, xcur)) - I;
    b = alpha*R*kron(xcur, xcur) - (1-alpha)*v;
    xn = A \ b;
    xn = max(xn,0);
    xn = xn ./ sum(xn);
    res = alpha*R*kron(xn, xn) + (1-alpha)*v - xn;
    if norm(res, 1) < tol
        break
    end
    xcur = xn;
end
if it == maxit
    warning('Maximum numer of iterations reached');
end