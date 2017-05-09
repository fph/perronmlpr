function [x, it, itn, xs, res] = optimistic(alpha, v, R, tol, maxit);
% new Perron vector method
% starts from x0 = 0

if not(exist('tol','var')) || isempty(eps)
    tol = sqrt(eps);
end
if not(exist('maxit','var')) || isempty(maxit)
    maxit = 10000;
end
n = length(v);
scaledR = alpha*R;
scaledv = (1-alpha)*v;

xs = [];
res = [];

[f, itn] = newton_from_zero(alpha, v, R, tol, maxit);
if sum(f) > 1 - tol
    warning 'is your problem supercritical?'
end
% beta = 1 - sum(f);
beta = 2- 1/ alpha;
y = beta * ones(n,1)/n;
for it = 1:maxit
    matrix = scaledR*(kron(f+y,eye(n)) + kron(eye(n),f));
%   [u, lambda] = perronvector(matrix, 'power',1.e-15,v); %lambda = 1 should hold
    %%%
    gamma=0.05;
      matrix1=gamma*matrix+(1-gamma)*eye(n);
  [U1,S1,V1]=svd(eye(n)-matrix1);
  u=V1(:,end);

      u = u / sum(u) * beta;
    % [v, lambda2] = perronvector(matrix', 'eig');  % v = 1 should hold
    e = ones(n,1);
%   gamma=0.99;
%      matrix=gamma*matrix+(1-gamma)*eye(n);
      
%%% con gamma=0.9999; alpha=0.99; R=R6_3, trova x =
%    -3.658595245222262e-02
%     2.641380990566147e-03
%     7.762161552420443e-03
%     1.133177544253783e+00
%    -4.478162428374238e-02
%    -6.221351006080490e-02
% con 5000 iterazioni

%%% con gamma=0.99999; alpha=0.99; R=R6_3, trova x =
%      1.999362874637594e-01 + 1.681924156479540e-04i
%      6.621960862221354e-03 + 1.511558792182352e-05i
%      1.164651232182554e-01 + 2.054978676157105e-04i
%      2.230890840471258e-01 - 7.613933911094400e-04i
%      7.997630648481438e-02 + 1.011120808573484e-04i
%      3.739112379238236e-01 + 2.714754390666035e-04i
% con 2100 iterazioni

%%% con gamma=0.9; alpha=0.99; R=R6_3, trova x =
%    -3.658595245223976e-02
%     2.641380990567006e-03
%     7.762161552425721e-03
%     1.133177544253843e+00
%    -4.478162428376340e-02
%    -6.221351006083253e-02
% con 46 iterazioni

% con gamma=0.99999, alpha=0.99; R=R6_3, usando matrix1 per la svd
% e matrix non rilassata per lo jacob, si trova x =
%    -3.208915248381875e-01
%     1.549892779870952e-01
%    3.123966330745115e-01
%     1.522991615201972e+00
%    -4.268877611090927e-01
%    -2.425982403162990e-01
% con 420 iterazioni

% con gamma=0.9999999, alpha=0.99; R=R6_3, usando matrix1 per la svd
% e matrix non rilassata per lo jacob, si trova x =
%
%     4.382072194627117e-02
%     2.224192630619705e-03
%     9.256490884021550e-03
%     8.191682635124664e-01
%     3.121744066976043e-02
%     9.431289035686066e-02
% con 7000 iterazioni
%%% stessa soluzione con gamma=0.1 e 2500 iter; gamma=0.05 e 980 iter

    %%%
%    Jac = (eye(n)- u*e'/(e'*u)) * pinv(matrix - eye(n)) * (u*e'/(e'*u) - eye(n)) * scaledR * kron(eye(n),u);
%    Jac = (eye(n)- u*e'/(e'*u)) * inv(eye(n)-matrix+u*e'/(e'*u)) * (scaledR*kron(eye(n),u)- alpha*u*e');
    Jac = inv(eye(n)-matrix+u*e'/(e'*u)) * (scaledR*kron(eye(n),u)) - alpha*u*e';
    [oldy, y] = deal(y, y - (eye(n) - Jac) \ (y - u));
    y = y / sum(y) * beta;
    if norm(oldy-y) < tol
        break
    end
    if nargout > 3
        xs = [xs f + y];
        x = y + f;
        res = [res norm(x - scaledR*kron(x,x) - scaledv)];
    end
end
x = y + f;
if it == maxit
    warning('Maximum numer of iterations reached');
end
