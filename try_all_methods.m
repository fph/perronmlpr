function [times, itns, validity] = try_all_methods(alpha, methods)
% Tries a set of methods on all the problems from Gleich's set
%
% [times, validity] = try_all_methods(alpha, methods)
%
% methods is a cell array of functions with signature result = f(alpha, v, R)
% alpha is the value of alpha to use.

experiments = load_tensor;

times = nan(length(experiments), length(methods));
itns = nan(length(experiments), length(methods));
validity = nan(length(experiments), length(methods));
for i = 1:length(experiments)
    experiment = experiments{i};
    R = load_tensor(experiment);
    v = 1/size(R,1) * ones(size(R,1), 1);
    for j = 1:length(methods)
        method = methods{j};
        try
            tic;
                [x, it] = method(alpha, v, R, sqrt(eps));
            times(i, j) = toc;
            itns(i, j) = it;
            validity(i, j) = check_solution(alpha, v, R, x, sqrt(eps));
        catch E
            times(i, j) = nan;
            validity(i, j) = nan;            
        end
    end
end