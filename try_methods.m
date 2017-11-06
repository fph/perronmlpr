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

[times, iters, results] = try_all_methods(0.99, methods);
sum(results==0)
% removes cases in which the algorithm produces a non-solution
iters_true = iters; iters_true(results~=0) = nan;
times_true = times; times_true(results~=0) = nan;

leg = {
    'FixedPoint', 
    'InnOut', 
    'Newton', 
    'Perron', 
%    'optimistic\_fixed', 
    'PerronNewton', 
%    'bootstrap\_newton', 
%    'bootstrap\_optimistic', 
%    'bootstrap\_optimistic\_newton', 
    'PerronNewton+Continuation',
    'Newton+Continuation',
    };

%[~, order] = sort(median(times,2));
order = 1:size(times,2);

set(0,'DefaultAxesFontSize',18);
set(0,'DefaultLegendFontSize',16);
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultLineMarkerSize', 10);

semilogy(times_true(order,:)); legend(leg)

legstring = join(leg,'\t'); legstring = legstring{1};
legstring = strrep(legstring, '\_', '_');
f = fopen('results.dat', 'wt');
fprintf(f, legstring);
fprintf(f, '\n');
dlmwrite('results.dat', times_true, 'delimiter', '\t', '-append');