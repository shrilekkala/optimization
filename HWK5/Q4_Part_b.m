n = 200;                % Dimension
x0 = 2 * ones(n, 1);    % Starting point of all 2's
tol = 1e-4;             % Tolerance for stopping criterion
max_iter = 1000;        % Maximum number of iterations

% Run the damped BFGS algorithm
[x_star, n_func_evals, grad_norms] = damped_BFGS(@rosenbrocknfgH, x0, tol, max_iter);

% Calculate empirical convergence rates
p = zeros(length(grad_norms)-2, 1);
for k = 2:length(grad_norms)-1
    p(k-1) = log(grad_norms(k+1) / grad_norms(k)) / log(grad_norms(k) / grad_norms(k-1));
end

% Compute the mean of p_k as an estimate for p
fprintf('An estimate for the rate of convergence is p = %.2f\n', mean(p));

% Create a semilogy plot
figure;
semilogy(1:length(grad_norms), grad_norms, 'o-');
grid on;
xlabel('Iteration');
ylabel('Norm of the gradient (log scale)');
title('Semilogy Plot of the Norm of Gradients');
