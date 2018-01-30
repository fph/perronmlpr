% This requires Bertini and BertiniLab

R = load_tensor('R6_3');
n = size(R, 1);

if not(exist('v', 'var'))
    v = 1/n * ones(n, 1);
end
    
names = "x" + [1:n];

polysyms(names{:});
eval('x=[' + names.join(';') + '];')

xx = [];
for i = 1:n
    xx = [xx; x(i)*x];
end

alphas = 0.7:0.01:0.8;

xs = [];

figure;

hold on;

for i = 1:length(alphas)
    alpha = alphas(i);

    % I tried to refactor this into functions, but scripts, functions and
    % BertiniLab gave trouble. :(
    
    % find solutions
    
    poly_system = BertiniLab('variable_group', x, 'function_def', alpha*R*xx + (1-alpha)*v - x);
    poly_system = poly_system.solve;
    sols = poly_system.match_solutions('real_finite_solutions');
    [x_re,~] = sols.x.real_imag;
    x_re = double(x_re);
    
    % filter solutions: some of them have sum=1, some have sum=
    
    
    plot(alpha*ones(size(x_re)), x_re, 'x');
end



hold off;