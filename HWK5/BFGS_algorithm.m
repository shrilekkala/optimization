function x = BFGS_algorithm(fgH, x0, tol, max_iter)
    % implements the BFGS algorithm 6.1 using backtracking line search
    % with Hessian approximation B_k (6.19)
    
    % fgH: Function that returns function value, gradient, and Hessian
    % x0: Starting point
    % tol: stopping criterion for gradient convergence
    % max_iter: maximum number of iterations for algorithm to run
    
    % Parameters for backtracking line search
    alpha_bar = 1;  % Initial step size
    rho = 0.25;      % Step size reduction factor
    c = 1e-4;       % Armijo rule constant

    % Initializations
    x = x0;
    iter = 0;
    I = eye(length(x0));    % Identity matrix
    B_k = I;                % Initial Hessian approximation B_0
    
    % Evaluate the function value and gradient at the starting point
    [~, g, ~] = fgH(x);
    
    % Main loop of the BFGS method
    while norm(g) > tol
        % Compute the search direction
        p_k = -B_k \ g;
        
        % Find the step length via backtracking
        alpha_k = backtracking_line_search(fgH, x, p_k, alpha_bar, rho, c);
        
        % Update the point
        x_new = x + alpha_k * p_k;
        
        % Calculate s_k and y_k for the Hessian update
        s_k = x_new - x;
        [~, g_new, ~] = fgH(x_new);
        y_k = g_new - g;
        
        % Update the Hessian approximation using (6.19)
        B_k = B_k - (B_k * (s_k * s_k') * B_k) / (s_k' * B_k * s_k) + (y_k * y_k') / (y_k' * s_k);

        
        % Update the iteration counter and the new points
        iter = iter + 1;
        x = x_new;
        g = g_new;

        % Check if maximum number of iterations have been reached
        if iter >= max_iter
            warning('Max number of iterations reached without convergence.')
            break
        end
    end
end