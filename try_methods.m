% Tests all methods to produce a plot of the CPU times.

relspeed = 0.01;
maxit=10000;

methods = {
    @fixed_point, 'FixedPoint';
    @innout, 'InnOut';
    @newton, 'Newton';
%    @(a,v,R,tol) tpr_method('innout',a,v,R,tol), 
%    @(a,v,R,tol) tpr_method('newton',a,v,R,tol),
    @optimistic, 'Optimistic';
%    @optimistic_fixed_point, 
    @optimistic_newton, 'Optimistic Newton';
    @(a,v,R,tol) bootstrap_first_derivative(@optimistic_newton,a,v,R,tol,maxit,relspeed), 'Bootstrap ON 1st derivative';
    @(a,v,R,tol) bootstrap_wrong_derivative(@optimistic_newton,a,v,R,tol,maxit,relspeed), 'Bootstrap ON wrong 1st';
    @(a,v,R,tol) bootstrap_second_derivative(@optimistic_newton,a,v,R,tol,maxit,relspeed), 'Bootstrap ON 2nd derivative';
    @(a,v,R,tol) bootstrap_derivativefree(@optimistic_newton,a,v,R,tol,maxit,relspeed), 'Bootstrap ON derivative-free';
%    @(a,v,R,tol) bootstrap_second_derivative_plus_previous(@optimistic_newton,a,v,R,tol,maxit,relspeed), 'Bootstrap ON derivative+prev';
    @(a,v,R,tol) bootstrap_first_derivative(@newton,a,v,R,tol,maxit,relspeed), 'Bootstrap N 1st derivative';
    @(a,v,R,tol) bootstrap_wrong_derivative(@newton,a,v,R,tol,maxit,relspeed), 'Bootstrap N wrong 1st';
    @(a,v,R,tol) bootstrap_second_derivative(@newton,a,v,R,tol,maxit,relspeed), 'Bootstrap N 2nd derivative';
    @(a,v,R,tol) bootstrap_derivativefree(@newton,a,v,R,tol,maxit,relspeed), 'Bootstrap N derivative-free';
    };

[times, iters, results] = try_all_methods(0.9, methods(:,1));
sum(results==0)
% removes cases in which the algorithm produces a non-solution
iters_true = iters; iters_true(results~=0) = nan;
times_true = times; times_true(results~=0) = nan;

set(0,'DefaultAxesFontSize',18);
set(0,'DefaultLegendFontSize',16);
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultLineMarkerSize', 10);
set(0,'defaultaxeslinestyleorder',{'--',':','-.'})

%semilogy(times_true(order,:)); legend(methods(:,2))

subplot(1,2,1);
title('iterations');
perfprof(iters_true); legend(methods(:,2), 'Location', 'southeast');
set(gca, 'XScale', 'log')
subplot(1,2,2);
title('CPU time');
perfprof(times_true); legend(methods(:,2), 'Location', 'southeast');
set(gca, 'XScale', 'log')

legstring = join(methods(:,2),'\t'); legstring = legstring{1};
legstring = strrep(legstring, ' ', '_');
legstring = strrep(legstring, '\_', '_');
f = fopen('results.dat', 'wt');
fprintf(f, legstring);
fprintf(f, '\n');
dlmwrite('results.dat', times_true, 'delimiter', '\t', '-append');