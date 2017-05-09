function [xn, it] = shifted(alpha, v, R, tol, maxit, x0)
% shifted iteration with gamma=1
% starts from x0 a stochastic vector
% Contains Gleich's code without the OO part

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end

if not(exist('x0','var')) || isempty(x0)
    x0 = (1-alpha)*v;
end

n = length(v);

gamma = 1; %contrary to the docstrings, the Gleich code solves with gamma=1 by default

Gamma = 1 / (1+gamma);
xcur = v;

for it=1:maxit
    % TODO make this iteration better
    y = alpha*(R*kron(xcur, xcur));
    z = y * Gamma + Gamma*(1-sum(y))*v;
    xn = z + (1-sum(z))*xcur;
    curres = alpha*R*kron(xn, xn) + (1-alpha)*v - xn;
    
    % switch solutions
    xcur = xn;
    
    % check for termination
    if norm(curres, 1) <= tol
        break;
    end
end

x = xcur ./ sum(xcur);
end
