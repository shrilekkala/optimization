function [x, n_func_evals, grad_norms] = damped_BFGS(fgH, x0, tol, max_iter)
    % Implements the Damped BFGS algorithm using backtracking line search
    
    % fgH: Function that returns function value, gradient, and Hessian
    % x0: Starting point
    % tol: Stopping criterion for gradient convergence
    % max_iter: Maximum number of iterations for algorithm to run

    % Parameters for backtracking line search
    alpha_bar = 1;  % Initial step size
    rho = 0.25;     % Step size reduction factor
    c = 1e-4;       % Armijo rule constant

    % Initializations
    x = x0;
    n_func_evals = 0;
    % Initialize a vector to 
    grad_norms = zeros(max_iter, 1); % Vector to store the norm of gradients
    iter = 0;
    I = eye(length(x0));    % Identity matrix
    B_k = I;                % Initial Hessian approximation B_0
    
    % Evaluate the function value and gradient at the starting point
    [~, g, ~] = fgH(x);

    % Update function evaluation count
    n_func_evals = n_func_evals + 1;
    
    % Main loop of the BFGS method
    while norm(g) > tol

        % Compute the search direction using B_k
        p_k = -B_k \ g;
        
        % Find the step length via backtracking
        [alpha_k, n_evals] = backtracking_line_search(fgH, x, p_k, alpha_bar, rho, c);
        % Update the point
        x_new = x + alpha_k * p_k;
        % Update function evaluation count
        n_func_evals = n_func_evals + n_evals;
        
        % Calculate s_k and y_k for the Hessian update
        s_k = x_new - x;
        [~, g_new, ~] = fgH(x_new);
        y_k = g_new - g;

        % Update function evaluation count
        n_func_evals = n_func_evals + 1;
        
        % Damping procedure
        % Compute theta_k
        if s_k' * y_k >= 0.2 * s_k' * B_k * s_k
            theta_k = 1;
        else
            theta_k = (0.8 * s_k' * B_k * s_k) / (s_k' * B_k * s_k - s_k' * y_k);
        end

        % Compute r_k
        r_k = theta_k * y_k + (1 - theta_k) * B_k * s_k;
        
        % Update the Hessian approximation B_k using (18.16)
        B_k = B_k - (B_k * (s_k * s_k') * B_k) / (s_k' * B_k * s_k) + (r_k * r_k') / (s_k' * r_k);
        
        % Update the iteration counter and the new points
        iter = iter + 1;
        x = x_new;
        g = g_new;

        % Update the list of norms
        grad_norms(iter) = norm(g);

        % Check if maximum number of iterations have been reached
        if iter >= max_iter
            warning('Max number of iterations reached without convergence.');
            return;
        end
    end

    % Report number of iterations for convergence
    fprintf('Damped BFGS Converged in %d iterations.\n', iter);

    % Update the list of grad norms as per n iterations
    grad_norms = grad_norms(1: iter);
end
