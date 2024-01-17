function [x, fval, iter, grad_norms, x_iters] = steepest_descent_backtracking(x0, f, grad_f, c)
    % Initialization
    alpha_bar = 1;        % Initial step length
    rho = 0.5;            % Contraction factor
    small_c = 1e-4;       % Sufficient decrease constant
    x = x0;               % Starting point
    iter = 0;             % Iteration counter
    x_iters = [];       % Array to store iterates
    grad_norms = [];      % Array to store gradient norms
    machine_precision = 1e-16; % Machine precision limit
    max_iter = 1000;      % Maximum number of iterations
    
    % Pre-calculate the convergence bound constant based on c
    conv_bound_const = (c - 1)^2 / (c + 1)^2;

    % Main loop
    while iter < max_iter
        iter = iter + 1;
        grad = grad_f(x);        % Compute gradient at current point
        grad_norms(iter) = norm(grad); % Store the gradient norm
        x_iters(:, iter) = x; % Store the current iterate
        p = -grad;               % Search direction is negative gradient (steepest descent)
        
        % Check stopping criterion based on Theorem 3.3
        if grad_norms(iter)^2 * conv_bound_const < machine_precision
            fval = f(x);         
            break;
        end
        
        % Backtracking line search
        alpha = alpha_bar;
        while f(x + alpha * p) > f(x) + small_c * alpha * grad' * p
            alpha = rho * alpha; % Contract the step size
        end
        
        % Update the point
        x = x + alpha * p;
        
    end
    
    if iter == max_iter
        fval = f(x);
    end
    
    % Output results
    fprintf('Steepest descent with backtracking line search completed.\n');
    fprintf('Iterations: %d\n', iter);
    fprintf('Optimal point: \n');
    disp(x);
end