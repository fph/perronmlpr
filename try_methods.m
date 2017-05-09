% Tests all methods to produce a plot of the CPU times.

methods = {
    @fixed_point,
    @innout,
    @newton,
%    @(a,v,R,tol) tpr_method('innout',a,v,R,tol), 
%    @(a,v,R,tol) tpr_method('newton',a,v,R,tol),
%    @newton, 
    @optimistic,
%    @optimistic_fixed_point, 
    @optimistic_newton,
%    @bootstrap_newton,
%    @bootstrap_optimistic,
%    @bootstrap_optimistic_newton,
    @bootstrap_optimistic_newton_second_order_stepsize,
    @bootstrap_newton_second_order_stepsize
    };

[times, iters, results] = try_all_methods(0.999, methods);
sum(results==0)
% removes cases in which the algorithm produces a non-solution
iters_true = iters; iters_true(results~=0) = nan;
times_true = times; times_true(results~=0) = nan;

leg = {
    'fixed', 
    'innout', 
    'newton', 
    'optimistic', 
%    'optimistic\_fixed', 
    'optimistic_newton', 
%    'bootstrap\_newton', 
%    'bootstrap\_optimistic', 
%    'bootstrap\_optimistic\_newton', 
    'bootstrap\_optimistic\_newton\_second\_order\_stepsize',
    'bootstrap\_newton\_second\_order\_stepsize',
    };

%[~, order] = sort(median(times,2));
order = 1:size(times,2);
semilogy(times_true(order,:)); legend(leg)
