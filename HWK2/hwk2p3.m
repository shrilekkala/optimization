function x = hwk2p3(x0, c)

    % Define the objective function and its gradient
    f = @(x) x(1)^2 + c * x(2)^2;
    grad_f = @(x) [2 * x(1); 2 * c * x(2)];

    % Call the steepest descent function with backtracking line search
    [x, fval, iter, grad_norms, x_iters] = steepest_descent_backtracking(x0, f, grad_f, c);
end