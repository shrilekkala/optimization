% Initializations
ns = 20:20:5000; % range of n's for the experiments
nmhm_systems_solved = [];
dogleg_systems_solved = [];
nmhm_times = [];
dogleg_times = [];

% Stopping Criteria
tol = 1e-3;     % tolerance for gradient norm convergence
max_iter = 100; % maximum number of iterations

% Measure total time of experiments when start
start = tic;

for n = ns
    fprintf('Running experiments for n = %d\n', n);
    x0 = 2 * ones(n, 1); % Initialize x to a vector of 2's

    % Time NMHM and record number of systems solved
    tic;
    [nmhm_result, nmhm_linear_solves] = NMHM(@rosenbrocknfgH, x0, tol, max_iter);
    nmhm_time = toc;

    % Time Dogleg and record number of systems solved
    tic;
    [dogleg_result, k, dogleg_linear_solves, dogleg_evals] = Trust_Region_Dogleg(@rosenbrocknfgH, x0, tol, max_iter);
    dogleg_time = toc;
    
    % Store times and numbers in lists
    nmhm_systems_solved(end+1) = nmhm_linear_solves;
    dogleg_systems_solved(end+1) = dogleg_linear_solves;
    nmhm_times(end+1) = nmhm_time;
    dogleg_times(end+1) = dogleg_time;

    % Stop experiments when total time elapsed is over 3 minutes
    if toc(start) > 180
        fprintf('Total script time exceeded 3 minutes, breaking loop n = %d\n', n);
        break;
    end
end

% Plot the results
figure;
computed_ns = ns(1:length(nmhm_times));

subplot(1, 2, 1);
plot(computed_ns, nmhm_systems_solved, 'b-', ...
    computed_ns, dogleg_systems_solved, 'r-','LineWidth', 1.5);
title('Number of Linear Systems Solved');
xlabel('Order n of Rosenbrock function');
ylabel('Number of systems solved');
legend('NMHM', 'Dogleg');
grid on;

subplot(1, 2, 2);
plot(computed_ns, nmhm_times, 'b-', ...
    computed_ns, dogleg_times, 'r-','LineWidth', 1.5);
title('Total Computation Time');
xlabel('Order n of Rosenbrock function');
ylabel('Compute time (s)');
legend('NMHM', 'Dogleg');
grid on;

