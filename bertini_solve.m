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

%alphas = 0.95:0.0002:1;
alphas = linspace(0.7,1,300);

xs = [];

figure;

hold on;


% copied from try_methods.m
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultLegendFontSize',16);
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultLineMarkerSize', 10);
set(0,'defaultaxescolororder',[0 0 0; 0 0 0.8; 1 0 0; 0 1 0]);
set(0,'defaultaxeslinestyleorder',{'-','--',':','-.'})

nsols = NaN;
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
    
    % filter solutions: some of them have sum=1, some have
    % sum=(1-alpha)/alpha
    
    x_re = x_re(:, sum(x_re) > mean([1, (1-alpha)/alpha]));
    x_re = x_re(:, all(x_re >= sqrt(eps)));
    
    x_re = x_re(1,:); % only plots first component
     
    if nsols ~= size(x_re, 2)
        fprintf('%d solutions for alpha=%g', size(x_re, 2), alpha);
    end
    nsols = size(x_re, 2);
    
    plot(alpha*ones(size(x_re)), x_re, 'xb');
end

xlabel('\alpha')

hold off;