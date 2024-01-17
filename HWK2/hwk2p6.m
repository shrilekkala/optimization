% Define the function f and its gradient for c = 1000
c = 1000;
f = @(x) x(1)^2 + c * x(2)^2;
grad_f = @(x) [2 * x(1); 2 * c * x(2)];

% Define starting points
starting_points = [1, 1; -1, 1; -1, -1; 1, -1];

% Initialize figure for plotting
figure;
hold on;

for i = 1:size(starting_points, 1)
    x0 = starting_points(i, :)';
    disp(['Starting Point: ', num2str(x0')])
    [x_opt, fval_opt, iter, grad_norms, x_iters] = steepest_descent_backtracking(x0, f, grad_f, c);
    
    % Calculate true error at each iteration
    % The minimum is at the origin, so we can just take the norm of each value
    true_errors = vecnorm(x_iters);
    
    % Plot the convergence metric (gradient norms)
    subplot(2, 4, i);
    hold on;
    semilogy(1:iter, grad_norms, 'DisplayName', ['x_0: ', num2str(x0')]);
    title('Convergence Metric ');
    xlabel('Iteration');
    ylabel('Log of Gradient Norm');
    
    % Plot the expected convergence from part 2
    conv_bound_const_sqrt = (c - 1) / (c + 1);
    semilogy(1:iter, grad_norms(1) * conv_bound_const_sqrt.^(0:iter-1),'-.','LineWidth', 1.5, 'DisplayName', 'Expected Rate');
    legend show;
    hold off;
    set(gca, 'yscale', 'log')
    
    % Plot the true error
    subplot(2, 4, i+4);
    semilogy(1:iter, true_errors, 'DisplayName', ['x_0: ', num2str(x0')]);
    title('True Error');
    xlabel('Iteration');
    ylabel('Log of True Error');
    legend show;
end

hold off;
