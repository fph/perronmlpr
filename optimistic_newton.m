function [x, it, itn, xs, ress] = optimistic_newton(alpha, v, R, tol, maxit, x0);
% Newton variant of the Perron method
%
% optional: solution guess x0

if not(exist('tol','var')) || isempty(tol)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end
n = length(v);
scaledR = alpha*R;
scaledv = (1-alpha)*v;

xs = [];
ress = [];

[f, itn] = newton_from_zero(alpha, v, R, tol, maxit);
if alpha == 0.5
    x = f;
    it = itn;
    xs = [x];
    res = alpha*R*kron(x, x) + (1-alpha)*v - x;
    return
end

if sum(f) > 1 - tol
    warning 'is your problem supercritical?'
end

beta = 2- 1/ alpha; % beta = 1 - sum(f);

if not(exist('x0','var')) || isempty(x0)
    y = beta * ones(n,1)/n;
else
    y = x0 - f;
    y = y / sum(y) * beta;
end

for it = 1:maxit-itn
    matrix = scaledR*(kron(f+y,eye(n)) + kron(eye(n),f));
    [u, lambda] = perronvector(matrix, 'eig'); %lambda = 1 should hold
    u = u / sum(u) * beta;
    % [v, lambda2] = perronvector(matrix', 'eig');  % v = 1 should hold
    e = ones(n,1);
%    Jac = (eye(n)- u*e'/(e'*u)) * pinv(matrix - eye(n)) * (u*e'/(e'*u) - eye(n)) * scaledR * kron(eye(n),u);
%    Jac = (eye(n)- u*e'/(e'*u)) * inv(eye(n)-matrix+u*e'/(e'*u)) * (scaledR*kron(eye(n),u)- alpha*u*e');
    Jac = (eye(n)-matrix+u*e'/(e'*u)) \ (scaledR*kron(eye(n),u)) - alpha*u*e';
   [oldy, y] = deal(y, y - (eye(n) - Jac) \ (y - u));
%   [oldy, y] = deal(y, y - (eye(n) - Jac + 0.01*eye(n)) \ (y - u)); %TODO testing
    y(y<0) = 0;
    y = y / sum(y) * beta;
    %%%%%%%%%%
%    w= [R*kron(eye(n),y);ones(1,n)]\[y;1] ; %%% se questo sistema sovradet ha soluz, lo Jac e' singolare
%    norm([R*kron(eye(n),y);ones(1,n)]*w-[y;1] )
    %%%%%%%%%%
    x = y + f; res = alpha*R*kron(x, x) + (1-alpha)*v - x;
    if nargout > 3
        xs = [xs x];
        ress = [ress res];
    end    
    if norm(res, 1) < tol
        break
    end
end
x = y + f;
it = itn + it;
if it == maxit
    warning('Maximum numer of iterations reached');
end
