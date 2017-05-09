function [x, it] = fixed_point(alpha, v, R, tol, maxit);
% fixed point iteration in "Multilinear Pagerank"
% starts from x0 = 0

if not(exist('eps','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

scaledv = (1-alpha)*v;
scaledR = alpha*R;
x = zeros(length(v),1);
for it = 1:maxit
    [oldx, x] = deal(x, scaledR*kron(x,x) + scaledv);
    if norm(oldx-x) < tol
        break
    end
end
