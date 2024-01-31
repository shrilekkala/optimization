function [x, n_linear_systems] = NMHM(fgH, x0, tol, max_iter)
    % Implements Newton's method with Hessian modification
    
    % fgH: Function that returns function value, gradient, and Hessian
    % x0: Starting point
    % tol: stopping criterion for gradient convergence
    % max_iter: maximum number of iterations for algorithm to run

    % Initialization
    x = x0;
    n_linear_systems = 0;
    
    % Parameters for backtracking line search
    alpha_bar = 1;  % Initial step size
    rho = 0.5;      % Step size reduction factor
    c = 1e-4;       % Armijo rule constant
    
    % Parameters for Hessian modification
    beta = 1e-3;    % Hessian modification factor
    
    % Algorithm 3.2 (Line Search Newton with Modification)
    for k = 1:max_iter

        % Calculate gradient and Hessian at current point
        [f_val, g, H] = fgH(x);
        
        % Modify the Hessian matrix using Algorithm 3.3
        Bk = modify_hessian_with_cholesky(H, beta);
        
        % Solve Bk * pk = - Grad(f(xk)) for pk
        pk = - Bk \ g;
        
        % Update linear system solve count
        n_linear_systems = n_linear_systems + 1;
        
        % Find step size using backtracking Algorithm 3.1
        alpha_k = backtracking_line_search(fgH, x, pk, alpha_bar, rho, c);
        
        % Update the point
        x = x + alpha_k * pk;
        
        % Check for convergence
        if norm(g) < tol
            fprintf('NMHM Converged in %d iterations.\n', k);
            return;
        end
    end
    warning('Max number of iterations reached without convergence.')
end