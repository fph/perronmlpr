function plot_minsol(R, v)
n = size(R, 1);

if not(exist('v', 'var'))
    v = 1/n*ones(n, 1);
end

sols = [];
alphas = linspace(0, 0.9, 50);
for i = 1:length(alphas)
    sols = [sols newton_from_zero(alphas(i), v, R)];
end
plot(alphas, sols);