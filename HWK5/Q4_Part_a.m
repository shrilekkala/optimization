% Parameters
tol = 1e-6;       % Stopping tolerance for gradient
max_iter = 1000;  % Maximum number of iterations

for n = 30:10:200
    fprintf('\nDimension: %d\n', n);

    % Starting point of all 2's
    x0 = 2 * ones(n, 1);
    
    % Run the damped BFGS algorithm
    [x_star, n_func_evals, ~] = damped_BFGS(@rosenbrocknfgH, x0, tol, max_iter);
    
    % Report the results
    fprintf('Number of Function Evaluations: %d\n', n_func_evals);
end
